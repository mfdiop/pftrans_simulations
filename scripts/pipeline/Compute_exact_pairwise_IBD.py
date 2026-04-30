# Compute exact pairwise IBD between all samples (by tree intervals)    
def pairwise_ibd(ts):
    samples = ts.samples()
    L = ts.sequence_length
    ibd = np.zeros((len(samples), len(samples)))

    for tree in ts.trees():
        interval = tree.interval  # length
        for i, a in enumerate(samples):
            for j, b in enumerate(samples):
                mrca = tree.mrca(a, b)
                # if they share the same ancestor in this tree, count the interval
                if mrca != tskit.NULL:
                    ibd[i,j] += interval

    ibd = ibd / L
    return ibd


# Corrected Version Here's a more accurate implementation:
import tskit
import numpy as np

def pairwise_ibd(ts, time_threshold=None, min_length=0):
    """
    Compute pairwise IBD segments between samples.
    
    Parameters:
    -----------
    ts : tskit.TreeSequence
        Input tree sequence
    time_threshold : float, optional
        Maximum TMRCA for considering segments as IBD (in generations)
        If None, uses all shared ancestry
    min_length : float
        Minimum segment length to count as IBD
    
    Returns:
    --------
    ibd : np.ndarray
        Matrix of IBD proportions (genome fraction shared)
    """
    samples = ts.samples()
    n = len(samples)
    L = ts.sequence_length
    ibd = np.zeros((n, n))
    
    for tree in ts.trees():
        span = tree.interval.right - tree.interval.left
        
        if span < min_length:
            continue
            
        for i, a in enumerate(samples):
            for j, b in enumerate(samples):
                if i == j:
                    continue  # Skip self-comparison
                
                mrca = tree.mrca(a, b)
                
                if mrca != tskit.NULL:
                    # Check time threshold if specified
                    if time_threshold is None:
                        ibd[i, j] += span
                    else:
                        mrca_time = tree.time(mrca)
                        if mrca_time <= time_threshold:
                            ibd[i, j] += span
    
    # Normalize by sequence length
    ibd = ibd / L
    return ibd