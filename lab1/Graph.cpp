#include "Graph.h"
#include "Net.h"
#include "Node.h"

Node *Graph::get_or_create_node(int index) {
    // Use node_map (std::map, O(log n)) for lookup instead of linear scan (O(n)).
    // For large graphs with hundreds of thousands of nodes, linear scan in a
    // growing vector would make this O(n^2), causing severe slowdown.
    auto it = node_map.find(index);
    if (it != node_map.end()) {
        return it->second;
    }
    Node *node = new Node(index);
    nodes.push_back(node);
    node_map[index] = node;
    return node;
}

Net *Graph::add_net(int index) {
    Net *net = new Net(index);
    this->nets.push_back(net);
    net_map[index] = net;
    return net;
}