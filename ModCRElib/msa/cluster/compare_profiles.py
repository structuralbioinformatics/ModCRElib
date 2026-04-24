import numpy as np

def padding(arr,final_length=50,mode="right"):
    """
    Pads or cuts a 2D numpy array (B, L) along the last dimension.

    Args:
        arr (np.ndarray): Input array of shape (B, L).
        final_length (int): Desired length after padding/cutting.
        mode (str): Padding mode: 'left', 'right', or 'center'.

    Returns:
        np.ndarray: Array of shape (B, final_length).
    """
    if arr.ndim != 2:
        raise ValueError("Input must be a 2D array of shape (B, L)")

    batch, length = arr.shape

    # Case 1: Cut if too long
    if length > final_length:
        if mode == "right":  # keep left part
            return arr[:, :final_length]
        elif mode == "left":  # keep right part
            return arr[:, -final_length:]
        elif mode == "center":  # cut equally from both sides
            start = (length - final_length) // 2
            end = start + final_length
            return arr[:, start:end]
        else:
            raise ValueError("mode must be 'left', 'right', or 'center'")

    # Case 2: Pad if too short
    pad_size = final_length - length
    if mode == "right":
        pad_left, pad_right = 0, pad_size
    elif mode == "left":
        pad_left, pad_right = pad_size, 0
    elif mode == "center":
        pad_left = pad_size // 2
        pad_right = pad_size - pad_left
    else:
        raise ValueError("mode must be 'left', 'right', or 'center'")

    return np.pad(arr, ((0, 0), (pad_left, pad_right)), mode="constant", constant_values=0)



def is_padding(column, threshold=0.0001):
    """
    Determine if a column corresponds to padding (all zeros or near-zero).
    """
    return column.sum() < threshold


def pearson_correlation(a, b,gap_penalty,threshold=0.0001):
    """
    Calculate Pearson correlation between two 1D-tensors.
    """
    if is_padding(a, threshold) or is_padding(b, threshold):
       return  -gap_penalty
    else:
       a_mean = a.mean()
       b_mean = b.mean()
       numerator =((a - a_mean) * (b - b_mean)).sum()
       denominator =np.sqrt(((a - a_mean) ** 2).sum()) * np.sqrt(((b - b_mean) ** 2).sum()) + 1e-8
       return numerator / denominator

def product(a,b,gap_penalty,threshold=0.0001):
    """
    Calculate dot product
    """

    if is_padding(a, threshold) or is_padding(b, threshold):
       return  -gap_penalty
    else:
       return (a * b).sum()

def kullback_leibler(a,b,gap_penalty,threshold=0.0001):
    """
    Calculate Kullback Leibler divergence
    """

    if is_padding(a, threshold) or is_padding(b, threshold):
       return  -gap_penalty
    else:
       return +(a * (np.log(a) - np.log(b))).sum()
    


def align_without_gaps(profile1, profile2, gap_penalty=0.5, threshold=0.0001):
    """
    Align two profiles without gaps by attempting all possible offsets.

    Args:
        profile1 (torch.Tensor): first profile of shape (4, L1).
        profile2 (torch.Tensor): Second profile of shape (4, L2).
        gap_penalty (float): Gap penalty for padding.
        threshold (float): minimum sum to be considered non-padding.

    Returns:
        max_score (float): The maximum alignment score.
        max_offset (int): The offset yielding maximum score.
        alignment (list of tuples): List of pairs (i, j) with matching positions.
    """
    L1 = profile1.shape[1]
    L2 = profile2.shape[1]

    max_score = float('-inf')
    max_offset = 0
    max_alignment = []

    # Try all possible offsets
    for offset in range(-L2+1, L1):
        score = 0
        alignment = []

        for i in range(L1):
            j = i - offset
            if j < 0 or j >= L2:
                continue

            col1 = profile1[:, i]
            col2 = profile2[:, j]

            #corr = pearson_correlation(col1, col2,gap_penalty,threshold)
            #corr = kullback_leibler(col1, col2,gap_penalty,threshold)
            corr = product(col1, col2,gap_penalty,threshold)
            score += corr
            alignment.append((i, j))

        if score > max_score:
            max_score = score
            max_offset = offset
            max_alignment = alignment.copy()

    return max_score, max_offset, max_alignment



# Example usage:

# Create two profiles with padding
profile1 = np.random.rand(4, 20)
profile2 = np.random.rand(4, 10)
profile3 = np.random.rand(4, 8)

profile1 = np.array(
 [[0.30017468 ,0.09420939 ,0.31260371 ,0.21689477 ,0.30879571 ,0.13937074 ,0.59566874 ,0.60783177 ,0.15212464 ,0.59759825 ,0.68761672 ,0.76494717 ,0.75355133 ,0.14071616 ,0.98248765 ,0.73724433 ,0.62294834 ,0.8732459  ,0.02084221 ,0.60096119],
 [0.92651669 ,0.57610969 ,0.32319033 ,0.64762615 ,0.26022946 ,0.60751543 ,0.11012558 ,0.98232935 ,0.5618884  ,0.59244851 ,0.86495986 ,0.58816366 ,0.87886669 ,0.01615445 ,0.15508351 ,0.04561064 ,0.38215413 ,0.52303797 ,0.76124283 ,0.75692127],
 [0.91946187 ,0.34880912 ,0.39552378 ,0.58855039 ,0.91610564 ,0.38277922 ,0.95529203 ,0.16562591 ,0.68729383 ,0.90817946 ,0.39706742 ,0.44126888 ,0.11225584 ,0.46808318 ,0.28530357 ,0.75088618 ,0.2656009  ,0.11094392 ,0.00795727 ,0.80912448],
 [0.73925381 ,0.01536006 ,0.52678562 ,0.91433777 ,0.75420477 ,0.73923654 ,0.41919009 ,0.71229939 ,0.85645588 ,0.20759987 ,0.51561869 ,0.17537583 ,0.17525748 ,0.77422442 ,0.13469913 ,0.73040476 ,0.6803141  ,0.8445494 ,0.16594266 ,0.4750031 ]]
 )

profile2 = np.array(
[[0.9051156  ,0.34969277 ,0.17959066 ,0.13525224 ,0.31515062 ,0.77092257 ,0.65733109 ,0.32749751 ,0.24253461 ,0.40252959],
 [0.44803286 ,0.48211524 ,0.52891759 ,0.81376805 ,0.62618829 ,0.55313191 ,0.79904126 ,0.25321574 ,0.2489375  ,0.41498628],
 [0.20980595 ,0.49111058 ,0.54203626 ,0.01672942 ,0.51838592 ,0.3739483 ,0.78670604 ,0.56553602 ,0.75790953 ,0.1905408 ],
 [0.47619524 ,0.96554684 ,0.40744212 ,0.9082903  ,0.51339929 ,0.63978421 ,0.65899225 ,0.43585935 ,0.20301949 ,0.84375403]]
)

profile3 = np.array(
[[0.33615287 ,0.38183845 ,0.76225093 ,0.7345508  ,0.50609033 ,0.33617003 ,0.14074239 ,0.94718676],
 [0.34038291 ,0.82427042 ,0.93195271 ,0.92249682 ,0.5565892  ,0.78072295 ,0.40391239 ,0.15677475],
 [0.31532874 ,0.93324838 ,0.64461175 ,0.08864189 ,0.30856656 ,0.76007311 ,0.08671389 ,0.41857227],
 [0.00562661 ,0.84290571 ,0.2302828  ,0.51226052 ,0.55438453 ,0.4070919 ,0.15673662 ,0.17215841]]
)

print("PROFILE1 ",profile1)
print("PROFILE2 ",profile2)
print("PROFILE3 ",profile3)
# Simulating padding by zeroing first and last columns

profile1=padding(profile1,final_length=50,mode="center")
#profile2=padding(profile1,final_length=50,mode="center")
#profile3=padding(profile1,final_length=50,mode="center")

#profile1[:, :5] = 0
#profile1[:, -5:] = 0

#profile2[:, :2] = 0
#profile2[:, -2:] = 0

#profile3[:, :8] = 0
#profile3[:, -8:] = 0

print("PROFILE1-padding ",profile1)
print("PROFILE2-padding ",profile2)
print("PROFILE3-padding ",profile3)

# Align
max_score, max_offset, max_alignment = align_without_gaps(profile1, profile2, gap_penalty=0.05)

print("Best alignment score =", max_score)
print("Best offset =", max_offset)
print("Alignment:")
print("Profile 1 -> Profile 2")
for i, j in max_alignment:
    print(f'{i} -> {j} score {product(profile1[:, i],profile2[:, j],gap_penalty=0.05)}')

max_score, max_offset, max_alignment = align_without_gaps(profile1, profile3, gap_penalty=0.05)

print("Best alignment score =", max_score)
print("Best offset =", max_offset)
print("Alignment:")
print("Profile 1 -> Profile 3")
for i, j in max_alignment:
    print(f'{i} -> {j} score {product(profile1[:, i],profile3[:, j],gap_penalty=0.05)}')

