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

def generate_lex_uc_with_events(n):
    """
    Run Hierholzer's algorithm but record a list of events to animate:
      - ('push', node)
      - ('take', from_node, to_node)
      - ('pop', node)
    After the traversal we append a sequence of ('highlight_cycle_edge', (u,v)) events
    that show the final cycle edges in order.
    Returns (uc_string, events, final_windows)
    """
    if n == 1:
        return "1", [('push', (1,))], [(1,)]
    if n == 2:
        return "12", [('push', (1,)), ('take', (1,), (2,)), ('push', (2,)), ('pop', (2,)), ('pop', (1,))], [(1,), (2,)]

    edges = list(itertools.permutations(range(1, n+1), n-1))
    adj = {}
    for edge in edges:
        v = edge[:-1]
        a = edge[-1]
        adj.setdefault(v, []).append(a)
    for v in adj:
        adj[v].sort()

    start = tuple(range(1, n-1))
    path = [start]
    edge_stack = []
    circuit = []
    events = []
    events.append(('push', start))

    # Hierholzer with events
    while path:
        v = path[-1]
        if adj[v]:
            a = adj[v].pop(0)
            u = v[1:] + (a,)
            edge_stack.append(a)
            events.append(('take', v, u))
            path.append(u)
            events.append(('push', u))
        else:
            popped = path.pop()
            events.append(('pop', popped))
            if edge_stack:
                circuit.append(edge_stack.pop())

    circuit.reverse()
    prefix = list(start)
    full_seq = prefix + circuit
    uc_str = ''.join(map(str, full_seq))

    # Build the windows (k = n-2) for the final cycle
    k = n - 2
    seq = [int(ch) for ch in uc_str]
    windows = []
    windows.append(tuple(seq[:k]))
    for i in range(len(seq) - k):
        win = windows[-1][1:] + (seq[i + k],)
        windows.append(win)

    # Build final cycle edges (pairs of consecutive windows)
    final_edges = []
    for i in range(len(windows) - 1):
        final_edges.append((windows[i], windows[i + 1]))
    # For completeness, the produced uc is of length n! so windows cover the full cycle.

    # append events to highlight the final cycle edges one-by-one
    for e in final_edges:
        events.append(('highlight_cycle_edge', e))

    return uc_str, events, windows

def display_shorthand_graph(n, uc_sequence="", interval=350):
    """
    Draws the directed graph of (n-2)-tuples → (n-2)-tuples.
    Animates according to the events produced by generate_lex_uc_with_events.
    """

    # Get uc and events
    uc, events, final_windows = generate_lex_uc_with_events(n)

    # Build vertex list (all (n-2)-tuples)
    k = n - 2
    vertices = list(itertools.permutations(range(1, n+1), k))
    N = len(vertices)

    # positions on unit circle
    pos = {
        v: (math.cos(2*math.pi*i / N), math.sin(2*math.pi*i / N))
        for i, v in enumerate(vertices)
    }

    # create plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_title(f"Shorthand-Permutation Graph for n={n}")
    ax.axis('off')

    # draw arrows for all possible edges (v -> u where u = v[1:]+(a,))
    edge_artists = {}
    for v, (x, y) in pos.items():
        missing = set(range(1, n+1)) - set(v)
        for a in missing:
            u = v[1:] + (a,)
            xu, yu = pos[u]
            arrow = ax.annotate(
                "",
                xy=(xu, yu), xytext=(x, y),
                arrowprops=dict(arrowstyle="->", lw=1, color="lightgray", shrinkA=5, shrinkB=5)
            )
            edge_artists[(v, u)] = arrow

    # nodes
    nodes = {}
    for i, (v, (x, y)) in enumerate(pos.items()):
        # initial coloring: white, except maybe the first vertex
        circ = ax.scatter(x, y, s=120, facecolor='white', edgecolor='black', zorder=3)
        ax.text(x, y + 0.08, str(v), ha='center', va='center', fontsize=8)
        nodes[v] = circ

    # Legend
    legend_patches = [
        mpatches.Patch(facecolor='green', edgecolor='black', label='On current path'),
        mpatches.Patch(facecolor='orange', edgecolor='black', label='Backtracked (popped)'),
        mpatches.Patch(facecolor='white', edgecolor='black', label='Unvisited')
    ]
    edge_patches = [
        mpatches.Patch(color='red', label='Edge taken'),
        mpatches.Patch(color='magenta', label='Final cycle edge')
    ]
    ax.legend(handles=legend_patches + edge_patches, loc='upper right', fontsize='small')

    def interpret_events_up_to(idx):
        """
        Simulate events[0:idx] while also tracking the edge-stack so that edges
        that are taken then later backtracked (popped) are returned as NOT taken.
        Returns:
        - current_path_nodes (set)
        - popped_nodes (set)
        - taken_edges (set)          # only edges currently active on the path
        - highlighted_cycle (list)   # final-cycle edges highlighted so far (preserve order)
        """
        path_stack = []
        popped = set()

        # taken_edges: set of edges currently 'active' (coloured red)
        taken_edges = set()
        # taken_stack: ordered stack of edges as they were taken, to simulate pops/backtrack
        taken_stack = []

        highlighted_cycle = []

        for e in events[:idx]:
            if e[0] == 'push':
                _, node = e
                path_stack.append(node)

            elif e[0] == 'take':
                # record the taken edge and push it onto the taken_stack
                _, a, b = e
                taken_stack.append((a, b))
                taken_edges.add((a, b))

            elif e[0] == 'pop':
                _, node = e
                # pop node from simulated path stack if it matches
                if path_stack and path_stack[-1] == node:
                    path_stack.pop()
                popped.add(node)

                # when a pop occurs, Hierholzer pops an edge from the edge stack too
                # if there is one; that edge should no longer be treated as "active".
                if taken_stack:
                    last_edge = taken_stack.pop()
                    if last_edge in taken_edges:
                        taken_edges.remove(last_edge)

            elif e[0] == 'highlight_cycle_edge':
                _, edge = e
                highlighted_cycle.append(edge)

        return set(path_stack), popped, taken_edges, highlighted_cycle


    # Animation update
    total_frames = len(events) + 10  # little hang at end
    def update(frame):
        # compute state up to this frame
        path_nodes, popped_nodes, taken_edges, highlighted_cycle = interpret_events_up_to(frame)

        # Reset visuals
        for v, circ in nodes.items():
            if v in path_nodes:
                circ.set_facecolor('green')
            elif v in popped_nodes:
                circ.set_facecolor('orange')
            else:
                circ.set_facecolor('white')

        # edges
        for (u, v), arrow in edge_artists.items():
            # default appearance
            arrow.arrow_patch.set_color('lightgray')
            arrow.arrow_patch.set_linewidth(1.0)

        # color edges that were taken as red
        for e in taken_edges:
            if e in edge_artists:
                edge_artists[e].arrow_patch.set_color('red')
                edge_artists[e].arrow_patch.set_linewidth(1.5)

        # highlight final cycle edges (magenta, thicker)
        for e in highlighted_cycle:
            if e in edge_artists:
                edge_artists[e].arrow_patch.set_color('magenta')
                edge_artists[e].arrow_patch.set_linewidth(2.5)
        # return artists (for blitting False is fine)
        artists = list(nodes.values()) + list(edge_artists.values())
        return artists

    ani = animation.FuncAnimation(
        fig, update, frames=total_frames,
        interval=interval, blit=False, repeat=False
    )
    plt.show()

if __name__ == "__main__":
    n = 4
    uc = generate_lex_uc(n)
    uc = uc[:math.factorial(n)]
    print("n=",n," →", uc)

    display_shorthand_graph(n, interval=700)
