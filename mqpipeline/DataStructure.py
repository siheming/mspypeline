from collections import defaultdict as ddict
from collections import deque
import pandas as pd


class DataNode:
    def __init__(self, name="", level=0, parent=None, data=None, children=None):
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

    def aggregate(self, method="mean", go_max_depth=False, index=None):
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


class DataTree:
    def __init__(self, root=None):
        self.root = root
        self.level_keys_full_name = ddict(list)
        # TODO self.level_keys_name = ddict(list)?
        # TODO self.methods ?

    def __getitem__(self, key, sep="_"):
        key_split = key.split(sep)
        start = self.root
        for k in key_split:
            start = start[k]
        return start

    @classmethod
    def from_analysis_design(cls, analysis_design, data=None, should_aggregate_technical_replicates=True):
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

    def add_data(self, data):
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