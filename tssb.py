import sys
import scipy.stats
import copy
from time import *
from numpy import *
from numpy.random import *
from util import *
from alleles import *
from copy import deepcopy, copy


class TSSB(object):
    min_dp_alpha = 0.00000000001
    max_dp_alpha = 1

    def __init__(self, dp_alpha, root_node=None, data=None, min_depth=0, max_depth=15):
        if root_node is None:
            raise Exception("Root node must be specified.")

        self.min_depth = 0
        self.max_depth = max_depth
        self.dp_alpha = dp_alpha
        self.data = data
        self.num_data = 0 if data is None else len(data)  # data.shape[0] #shankar
        self.root = root_node
        root_node.tssb = self

        self.assignments = []
        for n in range(self.num_data):
            self.root.add_datum(n)
            self.assignments.append(self.root)

    # Gibbs updates only
    # refer Equations (11) and (12) of treeCRP paper
    # Paper is unclear on what to do about internal nodes with only one datum.
    # Currently killing nodes that are singleton when their SSM comes up for resampling
    def resample_assignments(self):
        n_order = range(self.num_data)
        shuffle(n_order)
        for n in n_order:
            K = 0
            I = self.num_data
            # Remove datum
            self.assignments[n].remove_datum(n)
            self.cull()  # Kills previously singleton nodes
            (wts, nodes) = self.get_mixture()
            for node in nodes:
                # Add empty node to each existing one
                node.spawn()
                # Increment our counter of occupied nodes
                K += 1
            ps = []
            (wts, nodes) = self.get_mixture()
            for node in nodes:
                # Equation 11 or 12 applies depending on if the node has data assigned to it
                if node.num_data() > 0:
                    ps.append(log(node.num_data() / (self.dp_alpha + I - 1.0)) + node.logprob(self.data[n]))  # Eq11
                else:
                    ps.append(log(1 / (K + 1.0) * self.dp_alpha / (self.dp_alpha + I - 1.0)) + node.logprob(
                        self.data[n]))  # Eq12
            sum_ps = logsumexp(ps)
            ps = [exp(x - sum_ps) for x in ps]
            node_index = argmax(multinomial(1, ps))
            selected_node = nodes[node_index]
            selected_node.add_datum(n)
            self.assignments[n] = selected_node
            self.cull()

    # TODO
    # mixing weights need to be sampled from a Dirichlet (check...)
    # wts ~ Dirichlet(N_1+\alpha/K, N_2+\alpha/K, ..., N_K+\alpha/K)
    def get_mixture(self):
        def descend(root):
            nodes = [root]
            for i, child in enumerate(root.children()):
                child_nodes = descend(child)
                nodes.extend(child_nodes)
            return nodes

        nodes = descend(self.root)
        wts = []  # TODO...
        return (wts, nodes)

    def cull(self):
        def descend(node):
            for child in copy(node.children()):
                descend(child)
            if not node.has_data():
                node.kill()

        descend(self.root)

    def complete_data_log_likelihood(self):
        weights, nodes = self.get_mixture()
        llhs = []
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.data_log_likelihood())
        return sum(array(llhs))

    def complete_log_likelihood(self):
        weights, nodes = self.get_mixture()
        llhs = [self.dp_alpha_llh(self.dp_alpha)]
        for i, node in enumerate(nodes):
            if node.num_local_data():
                llhs.append(node.data_log_likelihood())
        return sum(array(llhs))

    def checkpi(self):
        spi = sum(self.getpi())
        if spi < 1 - 10e-6 or spi > 1 + 10e-6:
            raise Exception("PI INCONSISTENT. sum(pi) = %f" % spi)

    def getpi(self, tp=0, new=0):
        pi = []

        def descend(node, tp):
            for child in node.children():
                descend(child, tp)
            if new:
                p = node.pi1[tp]
            else:
                p = node.pi[tp]
            pi.append(p)

        descend(self.root, tp)
        pi = array(pi)
        pi.shape = len(pi),
        return pi

    def resample_alpha(self):
        upper = self.max_dp_alpha
        lower = self.min_dp_alpha
        llh = self.dp_alpha_llh(self.dp_alpha)
        llh_s = log(rand()) + llh
        while True:
            new_dp_alpha = (upper - lower) * rand() + lower
            new_llh = self.dp_alpha_llh(new_dp_alpha)
            if new_llh > llh_s:
                break
            elif new_dp_alpha < self.dp_alpha:
                lower = new_dp_alpha
            elif new_dp_alpha > self.dp_alpha:
                upper = new_dp_alpha
            else:
                break
        self.dp_alpha = new_dp_alpha

    def dp_alpha_llh(self, dp_alpha):
        llh = 0
        weights, nodes = self.get_mixture()
        for i, node in enumerate(nodes):
            contribution = betapdfln(float(len(node.data) + 0.0001) / (self.num_data + 0.0002), 1.0, dp_alpha)
            llh += contribution
        return llh

    def child_parent_split_merge(self):
        weights, nodes = self.get_mixture()
        N = len(nodes)
        node_index = argmax(multinomial(1, [1.0 / N] * N))
        if rand() > 0.5 or nodes[node_index] == self.root:
            # split
            hasting_ratio = (N + 1.0) / N
            source = nodes[node_index]
            source.spawn()
            dest = source._children[-1]
            for child in [x for x in source.children() if x != dest]:
                source.remove_child(child)
                dest.add_child(child)
                child._parent = dest
                dest.params = dest.params + child.params
            nids = list(source.data)
            if len(nids) == 0:
                return hasting_ratio
            nids.sort(key=lambda x: float(self.data[x].a[0]) / self.data[x].d[0], reverse=True)
            for ids in range(int(floor(len(nids) / 2.0))):
                source.remove_datum(nids[ids])
                dest.add_datum(nids[ids])
                self.assignments[nids[ids]] = dest

        else:
            # merge
            source = nodes[node_index]
            dest = source._parent
            data = copy(source.data)
            for ids in data:
                source.remove_datum(ids)
                dest.add_datum(ids)
                self.assignments[ids] = dest
            for child in copy(source._children):
                source.remove_child(child)
                dest.add_child(child)
                child._parent = dest
            source.kill()
            hasting_ratio = (N - 1.0) / N
        return hasting_ratio

    def leaf_sibling_split_merge(self):
        weights, nodes = self.get_mixture()
        if len(nodes) == 1:
            return log(0)  # Can't do sibling SM on root node
        nodes = [x for x in nodes if len(x.children()) == 0]
        N = len(nodes)
        node_index = argmax(multinomial(1, [1.0 / N] * N))
        Nsib = len(nodes[node_index]._parent.children()) - 1
        target_index = argmax(multinomial(1, [1.0 / (Nsib + 1)] * (Nsib + 1)))
        if target_index == 0:
            # split
            hasting_ratio = 2.0 * N * (Nsib + 1) / ((N + 1) * (Nsib + 1))
            source = nodes[node_index]
            source._parent.spawn()
            dest = source._parent.children()[-1]
            nids = list(source.data)
            if len(nids) == 0:
                return hasting_ratio
            nids.sort(key=lambda x: float(self.data[x].a[0]) / self.data[x].d[0], reverse=True)
            for ids in range(int(floor(len(nids) / 2.0))):
                source.remove_datum(nids[ids])
                dest.add_datum(nids[ids])
                self.assignments[nids[ids]] = dest
        else:
            # merge
            hasting_ratio = N * (Nsib + 1) / ((N - 1) * Nsib * 2.0)
            source = nodes[node_index]
            merge_targets = copy(source._parent.children())
            merge_targets.remove(source)
            dest = merge_targets[target_index - 1]
            data = copy(source.data)
            for ids in data:
                source.remove_datum(ids)
                dest.add_datum(ids)
                self.assignments[ids] = dest
            source.kill()
        return hasting_ratio

    def get_nodes(self):
        def descend(root):
            node = [root]
            for child in root.children():
                child_nodes = descend(child)
                node.extend(child_nodes)
            return node
        return descend(self.root)
