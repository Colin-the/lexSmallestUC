import collections
import itertools
import math
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.animation as animation

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
    
    ani = animation.FuncAnimation(
        fig, update, frames=len(windows),
        init_func=init, blit=True, interval=interval, repeat=False
    )
    plt.show()

if __name__ == "__main__":
    n = 4
    uc = generate_lex_uc(n)
    uc = uc[:math.factorial(n)]
    print("n=",n," →", uc)

    dot_output = display_shorthand_graph(n, uc_sequence=uc, interval=700)
    print(dot_output)
    
    # 123124132134214324314234
    # 123124132134214324314234

    # 123412351243124513241325134213452135214321452314231524153215342514352435142531453245312541354215432543154235412534152345
    # 123412351243124513241325134213452135214321452314231524153215342514352435142531453245312541354215432543154235412534152345




    