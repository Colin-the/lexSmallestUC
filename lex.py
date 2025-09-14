import collections
import itertools
import math
import random
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

def compute_force_layout(vertices, edge_list, radius=1.0, iterations=120, area_scale=1.0):
    """
    Very small, self-contained force-directed layout for small graphs.
    vertices: list of vertex keys
    edge_list: list of (u,v) edges (directed) to apply attractive forces along
    radius: target radius scale (final positions normalized to roughly this)
    iterations: number of relaxation steps (small graphs converge fast)
    """
    import random
    # initial positions on circle (good starting point)
    N = len(vertices)
    pos = {v: [math.cos(2*math.pi*i/N), math.sin(2*math.pi*i/N)] for i,v in enumerate(vertices)}

    # parameters (tweak for taste)
    k = math.sqrt(area_scale * (radius*radius) / max(1, N))  # ideal pairwise distance
    temperature = radius * 0.6

    for it in range(iterations):
        # repulsive forces
        disp = {v: [0.0, 0.0] for v in vertices}
        for i in range(N):
            vi = vertices[i]
            xi, yi = pos[vi]
            for j in range(i+1, N):
                vj = vertices[j]
                xj, yj = pos[vj]
                dx = xi - xj
                dy = yi - yj
                dist = math.hypot(dx, dy) + 1e-6
                # repulsive force magnitude
                fr = (k*k) / dist
                ux = dx / dist
                uy = dy / dist
                disp[vi][0] += ux * fr
                disp[vi][1] += uy * fr
                disp[vj][0] -= ux * fr
                disp[vj][1] -= uy * fr

        # attractive forces (along directed edges; treat as undirected for attraction)
        for (u, v) in edge_list:
            xu, yu = pos[u]
            xv, yv = pos[v]
            dx = xu - xv
            dy = yu - yv
            dist = math.hypot(dx, dy) + 1e-6
            fa = (dist*dist) / k
            ux = dx / dist
            uy = dy / dist
            # pull u towards v and v towards u (symmetric)
            disp[u][0] -= ux * fa
            disp[u][1] -= uy * fa
            disp[v][0] += ux * fa
            disp[v][1] += uy * fa

        # apply displacements with a cooling factor
        for v in vertices:
            dx, dy = disp[v]
            disp_len = math.hypot(dx, dy)
            if disp_len > 0:
                # limit by temperature
                scale = min(disp_len, temperature) / disp_len
                pos[v][0] += dx * scale
                pos[v][1] += dy * scale

        # cool
        temperature *= 0.95

    # normalize distances to roughly 'radius'
    # compute centroid and scale so average distance ~ radius
    cx = sum(p[0] for p in pos.values()) / N
    cy = sum(p[1] for p in pos.values()) / N
    avg_dist = sum(math.hypot(p[0]-cx, p[1]-cy) for p in pos.values()) / N
    if avg_dist > 1e-6:
        scale = radius / avg_dist
        for v in pos:
            pos[v][0] = (pos[v][0] - cx) * scale
            pos[v][1] = (pos[v][1] - cy) * scale

    return {v: tuple(pos[v]) for v in vertices}

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

def display_shorthand_graph(n, uc_sequence="", interval=350, fancy = False, filename = "traversal_animation.gif"):
    """
    Draws the directed graph of (n-2)-tuples → (n-2)-tuples.
    Animates according to the events produced by generate_lex_uc_with_events.
    """

    # Get uc and events
    uc, events, final_windows = generate_lex_uc_with_events(n)
    k = n - 2
    
    # create plot
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.set_title(f"Shorthand-Permutation Graph for n={n}")
    ax.axis('off')
    
    if n == 3 or fancy:  
        vertices = list(itertools.permutations(range(1, n+1), k))
        N = len(vertices)

        # choose layout method: for very small graphs use force-directed to spread nodes
        if N <= 6:
            # build edge list for attractive forces (treat missing elements same as earlier)
            tmp_edges = []
            for v in vertices:
                missing = set(range(1, n+1)) - set(v)
                for a in missing:
                    u = v[1:] + (a,)
                    tmp_edges.append((v, u))
            # small tweak: make radius slightly larger for n=3 so nodes are more spread
            base_radius = 1.8 if n == 3 else 1.4
            pos = compute_force_layout(vertices, tmp_edges, radius=base_radius, iterations=140)
        else:
            # default circular layout for larger graphs
            pos = { v: (math.cos(2*math.pi*i / N), math.sin(2*math.pi*i / N)) for i, v in enumerate(vertices) }

        # draw nodes
        nodes = {}
        for i, (v, (x, y)) in enumerate(pos.items()):
            circ = ax.scatter(x, y, s=120, facecolor='white', edgecolor='black', zorder=3)
            ax.text(x, y + 0.08, str(v), ha='center', va='center', fontsize=8)
            nodes[v] = circ

        # draw arrows but use curved connectionstyle to avoid exact overlaps
        from matplotlib.patches import FancyArrowPatch

        edge_artists = {}
        # group edges by pair to assign different radii when many edges share the same straight line
        pair_counts = {}
        for v in vertices:
            missing = set(range(1, n+1)) - set(v)
            for a in missing:
                u = v[1:] + (a,)
                pair_counts.setdefault((v,u), 0)
                pair_counts[(v,u)] += 1

        # build an ordering index per pair (so we can alternate curvature)
        pair_index = {}
        for (v, u) in pair_counts:
            pair_index[(v, u)] = 0

        for v, (x, y) in pos.items():
            missing = set(range(1, n+1)) - set(v)
            for a in missing:
                u = v[1:] + (a,)
                xu, yu = pos[u]

                # assign a curvature 'rad' depending on how many times this pair has been drawn
                idx = pair_index[(v, u)]
                pair_index[(v, u)] = idx + 1

                # alternate sign and scale by index to separate parallel edges (if any)
                sign = 1 if (idx % 2 == 0) else -1
                rad_base = 0.12  # adjust to taste
                rad = sign * (1 + idx//2) * rad_base

                arrow = FancyArrowPatch(
                    (x, y), (xu, yu),
                    connectionstyle=f"arc3,rad={rad}",
                    arrowstyle='-|>', mutation_scale=12,
                    linewidth=1.0,
                    color="lightgray",
                    zorder=1
                )
                ax.add_patch(arrow)
                edge_artists[(v, u)] = arrow

    else:
        vertices = list(itertools.permutations(range(1, n+1), k))
        N = len(vertices)

        # positions on unit circle
        pos = {
            v: (math.cos(2*math.pi*i / N), math.sin(2*math.pi*i / N))
            for i, v in enumerate(vertices)
        }

        

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
    ax.legend(handles=legend_patches + edge_patches, loc='upper left', bbox_to_anchor=(0.8, 1.12), bbox_transform=ax.transAxes, fontsize='small')

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

        # Reset visuals for nodes
        for v, circ in nodes.items():
            if v in path_nodes:
                circ.set_facecolor('green')
            elif v in popped_nodes:
                circ.set_facecolor('orange')
            else:
                circ.set_facecolor('white')

        # Helper to set color/linewidth for either FancyArrowPatch or annotate() result
        def set_arrow_look(arrow_obj, color, lw):
            # FancyArrowPatch has set_color / set_linewidth
            if hasattr(arrow_obj, "set_color") and hasattr(arrow_obj, "set_linewidth"):
                try:
                    arrow_obj.set_color(color)
                    arrow_obj.set_linewidth(lw)
                except Exception:
                    # defensive: some matplotlib builds behave slightly differently
                    pass
            # annotate() returns an Annotation with .arrow_patch attribute (a FancyArrowPatch)
            elif hasattr(arrow_obj, "arrow_patch") and arrow_obj.arrow_patch is not None:
                try:
                    arrow_obj.arrow_patch.set_color(color)
                    arrow_obj.arrow_patch.set_linewidth(lw)
                except Exception:
                    pass
            else:
                # last resort: try to set properties via setp (very defensive)
                try:
                    import matplotlib as mpl
                    mpl.pyplot.setp(arrow_obj, color=color, linewidth=lw)
                except Exception:
                    pass

        # Reset all edges to default appearance
        for (u, v), arrow in edge_artists.items():
            set_arrow_look(arrow, 'lightgray', 1.0)

        # Color edges that are currently taken (active on stack) as red
        for e in taken_edges:
            if e in edge_artists:
                set_arrow_look(edge_artists[e], 'red', 1.5)

        # Highlight final cycle edges (magenta, thicker)
        for e in highlighted_cycle:
            if e in edge_artists:
                set_arrow_look(edge_artists[e], 'magenta', 2.5)

        # return artists (nodes + edges). For annotations/fancy patches both are artists.
        artists = list(nodes.values()) + list(edge_artists.values())
        return artists


    ani = animation.FuncAnimation(
        fig, update, frames=total_frames,
        interval=interval, blit=False, repeat=False
    )
    fps = 1000.0 / interval  # frames per second (interval is in ms)

    try:
        # use PillowWriter
        writer = animation.PillowWriter(fps=fps)
        ani.save(filename, writer=writer, dpi=80)
    except Exception as e:
        print("Failed to save GIF:", e)


    plt.show()

if __name__ == "__main__":
    n = 6
    uc = generate_lex_uc(n)
    uc = uc[:math.factorial(n)]
    print("n=",n," →", uc)

    #display_shorthand_graph(n, interval=700)
