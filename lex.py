import collections
import itertools
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation
from itertools import permutations
from math import factorial

def generate_lex_uc(n):
    # Base cases for n=1 and n=2
    if n == 1:
        return "1"
    if n == 2:
        return "12"
    
    # Now for n > 2

    # 1) Build perfect-hash for all (n-1)-tuples
    edges = list(itertools.permutations(range(1, n+1), n-1))
    
    # 2) Build adjacency: vertex -> sorted list of outgoing symbols
    #    A vertex is a tuple of length (n-2)
    adj = {}
    for edge in edges:
        v = edge[:-1]  # first n-2 symbols
        a = edge[-1]   # last symbol
        adj.setdefault(v, []).append(a)
    # sort each list so that pop(0) gives the smallest
    for v in adj:
        adj[v].sort()

    # 3) Hierholzer’s with lexicographic choice
    start = tuple(range(1, n-1))
    path = [start]
    edge_stack = []    # symbols we took to get from each vertex in path
    circuit = []       # will collect the edges in reverse order
    
    while path:
        v = path[-1]
        if adj[v]:
            # take the smallest unused edge
            a = adj[v].pop(0)
            path.append(v[1:] + (a,))
            edge_stack.append(a)
        else:
            # backtrack: v has no more edges
            path.pop()
            if edge_stack:
                circuit.append(edge_stack.pop())
    
    # 4) circuit is edge‐labels in reverse; reverse it
    circuit.reverse()
    
    # 5) prepend the starting window
    prefix = list(start)
    full_seq = prefix + circuit
    return ''.join(map(str, full_seq))

def display_shorthand_graph(n, uc_sequence="", interval=300):
    """
    Draws the directed graph of (n-2)-tuples → (n-2)-tuples for shorthand permutations
    using pure matplotlib (no NetworkX).
    """

    if isinstance(uc_sequence, str):
        seq = [int(ch) for ch in uc_sequence]
    else:
        seq = list(uc_sequence)
    
    # Build the list of window-vertices visited:
    k = n - 2
    windows = []
    # first window: first k symbols
    windows.append(tuple(seq[:k]))
    for i in range(len(seq) - k):
        win = windows[-1][1:] + (seq[i + k],)
        windows.append(win)

        
    # 1) Generate vertices as all (n-2)-tuples
    vertices = list(itertools.permutations(range(1, n+1), n-2))
    N = len(vertices)
    
    # 2) Assign each vertex a position on the unit circle
    pos = {
        v: (math.cos(2*math.pi*i / N), math.sin(2*math.pi*i / N))
        for i, v in enumerate(vertices)
    }
    
    # 3) Create plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_title(f"Shorthand-Permutation Graph for n={n}")
    ax.axis('off')
    
    # 4) Draw edges as arrows
    edge_artists = {}
    for v, (x, y) in pos.items():
        missing = set(range(1, n+1)) - set(v)
        for a in missing:
            u = v[1:] + (a,)
            xu, yu = pos[u]
            arrow = ax.annotate(
                "",
                xy=(xu, yu), xytext=(x, y),
                arrowprops=dict(arrowstyle="->", lw=1, color="black", shrinkA=5, shrinkB=5)
            )
            edge_artists[(v, u)] = arrow
    
    # 5) Draw nodes and labels
    nodes = {}
    for v, (x, y) in pos.items():
        #print(x," ",y," ",v)
        if v == vertices[0]:
            circ = ax.scatter(x, y, s=100, color='green', edgecolor='black')
        else:
            circ =  ax.scatter(x, y, s=100, color='white', edgecolor='black')
        ax.text(x, y+0.075, str(v), ha='center', va='center', fontsize=8)
        nodes[v] = circ

    def init():
        for circ in nodes.values():
            circ.set_facecolor('white')
        for arrow in edge_artists.values():
            arrow.arrow_patch.set_color('black')
        return list(nodes.values()) + list(edge_artists.values())
    
    # Animation update: for each frame, highlight node and edge
    def update(frame):
        # Reset previous highlights
        if frame > 0:
            prev_node = windows[frame - 1]
            nodes[prev_node].set_facecolor('white')
            
            prev_edge = (windows[frame - 1], windows[frame])
            edge_artists[prev_edge].arrow_patch.set_color('black')
        
        # Highlight current node
        curr_node = windows[frame]
        nodes[curr_node].set_facecolor('green')
        
        # Highlight current edge (if not the first)
        if frame > 0:
            curr_edge = (windows[frame - 1], windows[frame])
            edge_artists[curr_edge].arrow_patch.set_color('red')
        
        return [nodes[curr_node], edge_artists.get((windows[frame - 1], windows[frame]), None)]
    
    # ani = animation.FuncAnimation(
    #     fig, update, frames=len(windows),
    #     init_func=init, blit=True, interval=interval, repeat=False
    # )
    plt.show()

def build_discovery_tree_and_cycle(n):
    """
    Returns:
      parent: dict mapping each (n-1)-tuple window u to its parent (n-2)-tuple v
      cycle:  list of symbols of the Eulerian circuit (length = n!)
    """
    symbols = list(range(1, n+1))
    k = n - 1
    vert_len = k - 1  # = n-2

    # 1) Build adjacency: each (n-2)-tuple -> sorted list of the missing symbols
    adj = {
        v: sorted([a for a in symbols if a not in v])
        for v in permutations(symbols, vert_len)
    }

    # 2) Initialize Hierholzer’s structures + parent‑tracking
    start = tuple(symbols[:vert_len])  # e.g. (1,2,...,n-2)
    stack      = [start]               # current path of vertices (each len n-2)
    edge_stack = []                    # symbols of outgoing edges we’ve taken
    circuit    = []                    # will collect edges (in reverse)
    seen_u     = set()                 # which (n-1)-windows we’ve discovered
    parent     = {}                    # parent[u] = v

    # 3) Traverse the graph
    while stack:
        v = stack[-1]
        if adj[v]:
            # take the lex‑smallest unused symbol from v
            a = adj[v].pop(0)
            u = v + (a,)           # this is a new (n-1)-tuple “window”
            if u not in seen_u:
                seen_u.add(u)
                parent[u] = v     # record discovery parent
            # step to the next vertex = last (n-2) entries of u
            stack.append(u[1:])
            edge_stack.append(a)
        else:
            # dead end: retreat
            stack.pop()
            if edge_stack:
                circuit.append(edge_stack.pop())

    # 4) Reverse to get the forward cycle of symbols
    cycle = circuit[::-1]
    assert len(cycle) == factorial(n), (
        f"Expected cycle length {factorial(n)}, got {len(cycle)}"
    )

    return parent, cycle


def plot_discovery_tree(parent):
    """
    Plots the discovery tree given a parent mapping child->parent of same window size.
    """
    # Find root (window with no parent)
    # For universal cycle, start window is the root
    root = next(w for w in parent.values() if w not in parent)

    # Compute depth of each node
    depth = {root: 0}
    print(depth)
    for child in parent:
        print("child:", child)
        depth[child] = depth[parent[child]] + 1

    # Group nodes by depth
    levels = {}
    for node, d in depth.items():
        levels.setdefault(d, []).append(node)
    for d in levels:
        levels[d].sort()

    # Assign positions
    positions = {}
    for d, nodes in levels.items():
        count = len(nodes)
        for i, node in enumerate(nodes):
            positions[node] = (i - (count-1)/2, -d)

    # Plot edges and nodes
    fig, ax = plt.subplots(figsize=(10, 6))
    for child, par in parent.items():
        x0, y0 = positions[par]
        x1, y1 = positions[child]
        ax.annotate("", xy=(x1, y1), xytext=(x0, y0),
                    arrowprops=dict(arrowstyle="->", lw=1))
    for node, (x, y) in positions.items():
        ax.scatter(x, y, s=100, color='black')
        ax.text(x, y, "".join(map(str, node)),
                ha='center', va='center', color='white', fontsize=8,
                bbox=dict(boxstyle="circle,pad=0.3", fc="black"))
    ax.axis('off')
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    n = 4
    uc = generate_lex_uc(n)
    uc = uc[:math.factorial(n)]
    print("n=",n," →", uc)

    #dot_output = display_shorthand_graph(n, uc_sequence=uc, interval=700)

    parent, cycle = build_discovery_tree_and_cycle(n)
    print(parent)
    plot_discovery_tree(parent)

    print(f"Lex‑smallest UC for (n−1)-perms of {{1..{n}}} has length {len(cycle)}")
    print(cycle)
    
    # 123124132134214324314234
    # 123 124 132 134 214 324 314 234

    # [3, 1, 2, 4, 1, 3, 2, 1, 3, 4, 2, 1, 4, 3, 2, 4, 3, 1, 4, 2, 3, 4, 1, 2]

    # 123124132134214324314234

    # 123412351243124513241325134213452135214321452314231524153215342514352435142531453245312541354215432543154235412534152345
    # 1234 1235 1243 1245 1324 1325 1342 1345 2135 2143 2145 2314 2315 2415 3215 3425 1435 2435 1425 3145 3245 3125 4135 4215 4325 4315 4235 4125 3415 2345


    


    # this is the out of order part of the tree
    # 2415 3215 3425 1435 2435 1425 3145 3245 3125 4135 4215 4325 4315 4235 4125 3415 2345



    # parent(2415) = e{2435, 3415}


    # 241 321 342 143 243 142 314 324 312 413 421 432 431 423 412 341 234
    # 124 132 234 143 243 142 143 243 123 134 142 243 143 234 124 134 234

    # parent(2415) = 2345 swap 4 and 1 then rase 1 to it's max (3)

    """
    Rule:
        1) if the string is in order then appily raise the right most non-max sysbol to it's max postion 
        2) if the string is not in order then decrement the left most out of order symbol

        "out of order" meaning a_i <= i + 1
    
        
    4125->3415->2345
    4125->2415->3415->1345->2345


    423 5

    4315->2435
    you rotate the first 3 symbols and then add 1 to the first one
    """
    