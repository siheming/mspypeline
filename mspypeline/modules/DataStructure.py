from collections import defaultdict as ddict
from collections import deque
import pandas as pd
from typing import Union, Callable, Dict, List, Sized, Optional


class DataNode:
    def __init__(self, name: str = "",
                 level: int = 0,
                 parent: "DataNode" = None,
                 data: pd.Series = None,
                 children: Dict[str, "DataNode"] = None):
        """
        Default parameters will return a root node

        Parameters
        ----------
        name
            Name of the node
        level
            depth of the node
        parent
            Parent of this node
        data
            Is None when there are nodes below this one, which were not aggregated as technical replicated
        children
            Maps name of a child to a child node

        See Also
        --------
        DataTree: A class to help construct a node structure from data

        """
        self.name = name
        self.parent = parent
        self.data = data
        self.level = level
        self.children = {} if children is None else children
        if self.parent is not None:
            prefix = self.parent.full_name + "_" if self.parent.full_name != "" else ""
            self.full_name = prefix + self.name
        else:
            self.full_name = ""
        self.iter_counter = -1

    def __str__(self):
        return f"level {self.level}, name: {self.full_name}, n children: {self.get_total_number_children()}"

    def __repr__(self):
        return f"DataNode(name={self.name}, level={self.level}, parent={self.parent}, data={self.data}, children={self.children}"

    def __getitem__(self, key):
        return self.children[key]

    def __setitem__(self, key, node):
        self.children[key] = node

    def __iter__(self):
        self.iter_counter = -1
        return self

    def __next__(self):
        self.iter_counter += 1
        if self.iter_counter < len(self.children.keys()):
            idx = list(self.children.keys())[self.iter_counter]
            return self.children[idx]
        raise StopIteration

    def get_total_number_children(self, go_max_depth: bool = False) -> int:
        """
        Gets the number of children containing data below this node. If go_max_depth will search for the deepest
        DataNodes.

        Parameters
        ----------
        go_max_depth
            default false

        Returns
        -------
        int
            The number of all children below this node

        """
        n_children = 0
        if self.data is not None and not go_max_depth:
            return n_children
        queue = deque([self])
        while queue:
            parent = queue.popleft()
            for child in parent:
                if child.data is not None and child.children and go_max_depth:
                    queue += [child]
                elif child.data is None:
                    queue += [child]
                else:
                    n_children += 1
        return n_children

    def aggregate(self,
                  method: Union[None, str, Callable] = "mean",
                  go_max_depth: bool = False,
                  index: Optional[Union[str, pd.Index]] = None):
        """

        Parameters
        ----------
        method
            If None no aggregation will be applied. Otherwise needs to be accepted by pd.aggregate.
        go_max_depth
            If technical replicates were aggregated, this can be specified to use the unaggregated values instead.
        index
            Index to subset the data with. If None no index is applied

        Returns
        -------
        Union[pd.Series, pd.DataFrame]
            Result of the aggregation
        """
        queue = deque([self])
        data = []
        while queue:
            parent = queue.popleft()
            should_go_deeper = go_max_depth and parent.children
            if parent.data is not None and not should_go_deeper:
                if index is not None:
                    # append only the items in the index
                    series_data = parent.data.loc[index]
                    if not isinstance(series_data, pd.Series):
                        series_data = pd.Series(series_data, name=parent.data.name, index=[index])
                    data.append(series_data)
                else:
                    data.append(parent.data)
            else:
                for child in parent:
                    queue += [child]
        data = pd.concat(data, axis=1)
        if method is not None:
            data = data.aggregate(method, axis=1).rename(self.full_name, axis=1)
        return data

    def groupby(
            self, method: Union[str, Callable] = "mean",
            go_max_depth: bool = False,
            index: Union[None, str, pd.Index] = None
    ) -> Union[pd.Series, pd.DataFrame]:
        """
        consider each child a group then aggregate all children

        Parameters
        ----------
        method
            Will be passed to aggregate.
        go_max_depth
            Will be passed to aggregate.
        index
            Will be passed to aggregate.

        Returns
        -------
        data
            Result of the grouping

        See Also
        --------
        aggregate : Will be called on each of the groups

        """
        data = {child.name: child.aggregate(method, go_max_depth, index) for child in self}
        data = pd.concat(data, axis=1)
        new_col_names = [self.full_name]
        if method is None:
            new_col_names.append("level_1")
        data.columns = data.columns.set_names(new_col_names)
        return data


class DataTree:
    """Summary line

    Attributes
    ----------
    root: DataNode
        Does stuff
    level_keys_full_name: Dict[int, List[str]]
        Has all DataNode.full_name of a depth level of all levels

    """
    def __init__(self, root: DataNode):
        """

        Parameters
        ----------
        root
            The root node of the Tree.
        """
        self.root = root
        self.level_keys_full_name = ddict(list)
        # TODO self.level_keys_name = ddict(list)?
        # TODO self.methods ?

    def __getitem__(self, key: str, sep: str = "_"):
        # TODO maybe this should be moved to the data node
        key_split = key.split(sep)
        start = self.root
        for k in key_split:
            start = start[k]
        return start

    @classmethod
    def from_analysis_design(cls,
                             analysis_design: dict,
                             data: Union[None, pd.DataFrame] = None,
                             should_aggregate_technical_replicates: bool = True):
        """

        Parameters
        ----------
        analysis_design
            nested dict
        data
            Will be passed to add_data. If None no data is added
        should_aggregate_technical_replicates
            If True the lowest level of the analysis design is considered as a technical replicate and averaged

        Returns
        -------
        cls

        See Also
        --------
        add_data: will be called if data is not None
        aggregate_technical_replicates: will be called if should_aggregate_technical_replicates

        """
        root = DataNode()
        c = cls(root)
        queue = deque([(0, root, analysis_design)])
        while queue:
            level, parent, d = queue.popleft()
            if isinstance(d, dict):
                for child in d:
                    node = DataNode(name=child, level=level, parent=parent)
                    parent[child] = node
                    queue += [(level + 1, node, d[child])]
                    c.level_keys_full_name[level].append(node.full_name)
                    # c.level_keys_name[level].append(node.name)
        if data is not None:
            c.add_data(data)
        if should_aggregate_technical_replicates:
            c.aggregate_technical_replicates()
        return c

    def aggregate(
            self, key: Optional[str] = None,
            method: Union[None, str, Callable] = "mean",
            go_max_depth: bool = False,
            index: Optional = None
    ) -> Union[pd.Series, pd.DataFrame]:
        """

        Parameters
        ----------
        key
        method
        go_max_depth
        index

        Returns
        -------

        """
        if key is None:
            return self.root.aggregate(method, go_max_depth, index)
        elif isinstance(key, str):
            return self[key].aggregate(method, go_max_depth, index)
        else:
            raise ValueError(f"Invalid input for key: {key}, with type: {type(key)}")

    def groupby(self,
                key_or_index: Union[None, str, int] = None,
                new_col_name: str = None,
                method: Union[None, str, Callable] = "mean",
                go_max_depth: bool = False,
                index=None) -> Union[pd.Series, pd.DataFrame]:
        """

        Parameters
        ----------
        key_or_index
        new_col_name
        method
        go_max_depth
        index

        Returns
        -------

        """
        if key_or_index is None:
            return self.root.groupby(method, go_max_depth, index)
        elif isinstance(key_or_index, str):
            return self.root[key_or_index].groupby(method, go_max_depth, index)
        elif isinstance(key_or_index, int):
            data = {
                child_name: self[child_name].aggregate(method, go_max_depth, index)
                for child_name in self.level_keys_full_name[key_or_index]
            }
            data = pd.concat(data, axis=1)
            if new_col_name is None:
                new_col_names = ["level_0"]
            else:
                new_col_names = [new_col_name]
            if method is None:
                new_col_names.append("level_1")
            data.columns = data.columns.set_names(new_col_names)
            return data
        else:
            raise ValueError(f"Invalid input for key: {key_or_index}, with type: {type(key_or_index)}")

    def add_data(self, data: pd.DataFrame):
        """

        Parameters
        ----------
        data
            Data which will be used to fill the nodes with a Series. The column names of the data need to be the same as
            the full names of the DataNode.

        """
        queue = deque([self.root])
        while queue:
            parent = queue.popleft()
            for child in parent:
                queue += [child]
                if child.full_name in data.columns:
                    child.data = data.loc[:, child.full_name]

    def aggregate_technical_replicates(self):
        """
        Aggregates the deepest level to one level above by using aggregate

        """
        queue = deque([self.root])
        while queue:
            parent = queue.popleft()
            for child in parent:
                if child.children:
                    queue += [child]
                else:
                    parent.data = parent.aggregate()
