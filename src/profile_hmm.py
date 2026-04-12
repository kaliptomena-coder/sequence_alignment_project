def generate_profile(msa):
    """Building a position-specific frequency matrix from an existing MSA."""
    length = len(msa[0])
    profile = []
    for col in range(length):
        counts = {}
        for seq in msa:
            char = seq[col]
            # Counting occurrences of each amino acid or gap
            counts[char] = counts.get(char, 0) + 1

        # Converting counts to percentages
        freqs = {k: round(v / len(msa), 2) for k, v in counts.items()}
        profile.append(freqs)
    return profile

if __name__ == "__main__":
    # Using your Refined MSA result
    msa_result = [
        "-VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFP-T",
        "-VLSAADKTNVKAAWSKVGGHAGEYGAEALERMFLGFP-T",
        "-VHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPW-T"
    ]

    print("Profile HMM: Column Frequency Analysis")
    profile = generate_profile(msa_result)

    # Displaying the first 10 columns of the profile
    print(f"{'Column':<5} | {'Amino Acid Frequencies'}")
    print("-" * 35)
    for i in range(10):
        print(f"{i:<5} | {profile[i]}")