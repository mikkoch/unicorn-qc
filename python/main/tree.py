from main.gmm import *
from typing import *


class Node:
    def __init__(
            self,
            name:str,
            n_populations:int,
            n_samples: int,
            samples: List[str],
            is_leaf: bool,
            gmm: GMM=None):
        self.name = name
        self.n_populations = n_populations
        self.n_samples = n_samples
        self.samples = samples
        self.gmm = gmm
        self.is_leaf = is_leaf

    @staticmethod
    def is_root(
        branch:str,
        root_label:str='root') -> bool:
        return branch == root_label

    @staticmethod
    def get_parent(
        branch:str,
        separator:str='-')->str:
        parent_branch = separator.join(branch.split(separator)[:-1])
        if branch == parent_branch:
            raise Warning(
                "The parent of current branch {branch} is itself".format(branch=branch))
        return parent_branch

    @staticmethod
    def get_ancestors(
        branch:str,
        separator:str='-'
    )->List[str]:
        ancestors = list()
        while Node.get_parent(branch, separator) != branch:
            parent = Node.get_parent(branch=branch, separator=separator)
            branch = parent
            ancestors.append(parent)
        return ancestors

    def get_name(self)->str:
        return self.name

    def get_classifier(self)->GMM:
        return self.gmm

    def get_n_populations(self)->int:
        return self.n_populations

    def get_n_samples(self)->int:
        return self.n_samples

    def get_samples(self)->List[str]:
        return self.samples



class Tree:
    def __init__(self):
        self.__children__ = dict()

    def add_node(self, n: Node):
        node_name = n.get_name()
        children = self.__children__
        children.update({node_name: n})
        self.__children__ = children

    def get_node(self, name: str)->Node:
        return self.__children__.get(name)

    def get_all_nodes(self)->Dict[str, Node]:
        return self.__children__

    def get_all_node_names(self) -> List[str]:
        return list(self.__children__.keys())

    def get_leaves(self) -> Dict[str, Node]:
        return dict([(name, node) for name, node in self.__children__.items() if node.is_leaf])

    def get_assignment(
            self,
            output_sample_col:str='s',
            output_col:str='cluster_assignment',
            n_samples_thresh:int=100,
            hail_output:bool=True) -> Union[hl.Table, pd.DataFrame]:

        assign_pd_frames = list()
        for name, node in self.get_leaves().items():
            if node.get_n_samples() < n_samples_thresh:
                continue
            assign_pd = pd.DataFrame({
                output_sample_col:node.get_samples(),
                output_col: name})
            assign_pd_frames.append(assign_pd)

        out_assign_pd = pd.concat(assign_pd_frames)
        print(out_assign_pd)
        if hail_output:
            out_assign_ht = hl.Table.from_pandas(out_assign_pd, key=output_sample_col)
            return out_assign_ht
        else:
            return out_assign_pd
