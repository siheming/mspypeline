from collections import defaultdict as ddict
from collections import deque
import pandas as pd
from typing import Union, Callable, Dict


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
        level
        parent
        data
        children
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

    def __repr__(self):
        return f"level {self.level}, name: {self.full_name}, n children: {len(self.children)}"

    def __getitem__(self, key):
        return self.children[key]

    def __setitem__(self, key, node):
        self.children[key] = node

    def __iter__(self):
        return (node for node in self.children.values())

    def __next__(self):
        for child in self.__iter__():
            return child

    def aggregate(self,
                  method: Union[None, str, Callable] = "mean",
                  go_max_depth: bool = False,
                  index=None):
        """

        Parameters
        ----------
        method
            If None no aggregation will be applied. Otherwise needs to be accepted by pd.aggregate.
        go_max_depth
            If technical replicates were aggregated, this can be specified to use the unaggregated values instead.
        index
            Index to subset the data with.

        Returns
        -------
        Union[pd.Series, pd.DataFrame]
            Result of the aggregation
        """
        queue = deque([self])
        data = []
        while queue:
            parent = queue.popleft()
            has_child = next(parent) is not None
            should_go_deeper = go_max_depth and has_child
            if parent.data is not None and not should_go_deeper:
                # TODO maybe use the index here already
                data.append(parent.data)
            else:
                for child in parent:
                    queue += [child]
        data = pd.concat(data, axis=1)
        if index is not None:
            data = data.loc[index, :]
        if method is not None:
            data = data.aggregate(method, axis=1).rename(self.full_name, axis=1)
        return data

    # TODO
    def groupby(self, method="mean", go_max_depth=False, index=None):
        """
        consider each child a group then aggregate all children
        Parameters
        ----------
        method: Union[str, Callable]
            Will be passed to aggregate.
        go_max_depth: bool
            Will be passed to aggregate.
        index:
            Will be passed to aggregate.

        Returns
        -------
        Union[pd.Series, pd.DataFrame]
            Result of the grouping

        See Also
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
    level_keys_full_name: ddict(list)
        Has all DataNode.full_name of a depth level of all levels

    """
    def __init__(self, root):
        """

        Parameters
        ----------
        root: DataNode
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
                             analysis_design,
                             data: pd.DataFrame = None,
                             should_aggregate_technical_replicates: bool = True):
        """

        Parameters
        ----------
        analysis_design
        data
            Will be passed to add_data
        should_aggregate_technical_replicates

        Returns
        -------
        cls

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

    def aggregate(self,
                  key: Union[None, str] = None,
                  method: Union[None, str, Callable] = "mean",
                  go_max_depth: bool = False,
                  index=None):
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
                index=None):
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

    def add_data(self, data):
        """

        Parameters
        ----------
        data : pd.DataFrame
            Data which will be used to fill the nodes with a Series.

        """
        queue = deque([self.root])
        while queue:
            parent = queue.popleft()
            for child in parent:
                queue += [child]
                if child.full_name in data.columns:
                    child.data = data.loc[:, child.full_name]

    def aggregate_technical_replicates(self):
        queue = deque([self.root])
        while queue:
            parent = queue.popleft()
            if parent.children:
                if next(parent).children:
                    for child in parent:
                        queue += [child]
                else:
                    parent.data = parent.aggregate()
