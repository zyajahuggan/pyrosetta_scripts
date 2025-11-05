def identify_secondary_structure_spans(ss):
    blocks = []
    current = None
    start = None

    for i, j in enumerate(ss, start=1):
        # If we are entering a structured segment
        if j in ('E', 'H'):
            if current is None:
                current = j
                start = i
            elif j != current:
                blocks.append((start, i-1))
                current = j
                start = i
        else:  # j == 'L'
            if current is not None:
                blocks.append((start, i-1))
                current = None
                start = None
    # Close final block if necessary
    if current is not None:
        blocks.append((start, len(ss)))
    return blocks

def fold_tree_from_dssp_string(ss_string):
    N = len(ss_string)

    # --- Identify SSE blocks (runs of E/H) ---
    sse = []
    current = None
    start = None
    for i, ch in enumerate(ss_string, start=1):
        if ch in ('E','H'):
            if current is None:
                current = ch
                start = i
            elif ch != current:
                sse.append((start, i-1))
                current = ch
                start = i
        else:
            if current is not None:
                sse.append((start, i-1))
                current = None
                start = None
    if current is not None:
        sse.append((start, N))

    if not sse:
        raise ValueError("No helices or strands detected.")

    # --- Compute midpoints ---
    def mid(a,b): return (a+b)//2
    sse_mid = [mid(a,b) for (a,b) in sse]

    # --- Identify loops between SSEs ---
    loops = []
    loop_mid = []
    for (a1,b1),(a2,b2) in zip(sse, sse[1:]):
        if a2 - b1 > 1:
            loops.append((b1+1, a2-1))
            loop_mid.append(mid(b1+1, a2-1))

    root = sse_mid[0]  # anchor midpoint

    ft = FoldTree()
    jump_id = 1

    # --- First SSE: back to 1, forward to SSE1 end ---
    a0,b0 = sse[0]
    if root > 1:
        ft.add_edge(root, 1, -1)
    if root < b0:
        ft.add_edge(root, b0, -1)

    # --- Interleave each loop and the next SSE ---
    rows = max(len(loops), len(sse)-1)
    for k in range(rows):

        # Loop branch
        if k < len(loops):
            lm = loop_mid[k]; la,lb = loops[k]
            ft.add_edge(root, lm, jump_id); jump_id += 1
            if lm > la: ft.add_edge(lm, la, -1)
            if lm < lb: ft.add_edge(lm, lb, -1)

        # Next SSE branch
        if (k+1) < len(sse):
            sm = sse_mid[k+1]; sa,sb = sse[k+1]
            ft.add_edge(root, sm, jump_id); jump_id += 1
            forward_end = N if (k+1)==len(sse)-1 else sb
            if sm > sa: ft.add_edge(sm, sa, -1)
            if sm < forward_end: ft.add_edge(sm, forward_end, -1)

    try:
        ft.reorder(root)
    except Exception:
        pass

    return ft

def fold_tree_from_ss(pose):
    ss = Dssp(pose).get_dssp_secstruct()
    ss = ss[:pose.total_residue()].ljust(pose.total_residue(), 'L')
    return fold_tree_from_dssp_string(ss)

