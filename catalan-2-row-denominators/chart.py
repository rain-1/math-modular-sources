import csv
import matplotlib.pyplot as plt

def load_csv(filename):
    """Utility to load a CSV into a dictionary of column lists."""
    data = {}
    try:
        with open(filename, 'r') as f:
            reader = csv.reader(f)
            headers = next(reader)
            for h in headers:
                data[h] = []
            for row in reader:
                for h, val in zip(headers, row):
                    try:
                        val = float(val)
                    except ValueError:
                        if val == 'True': val = True
                        elif val == 'False': val = False
                    data[h].append(val)
    except FileNotFoundError:
        print(f"Warning: {filename} not found. Ensure it is in the same directory.")
    return data

def main():
    # Load the data from CSVs
    denom_data = load_csv("zudilin_denominators.csv")
    zudilin_data = load_csv("zudilin_rows.csv")
    nesterenko_data = load_csv("nesterenko_rows.csv")
    combined_data = load_csv("combined_rows.csv")

    # =========================================================================
    # PLOT 1: p=2 Denominator Valuation Bound
    # =========================================================================
    if denom_data:
        plt.figure(figsize=(7, 5))
        m = denom_data['m']
        v2 = denom_data['v2_denom']
        e_m = denom_data['e_m_bound']
        
        plt.plot(m, e_m, 'r--', label='Theoretical Bound ($e_m$)')
        plt.plot(m, v2, 'bo-', label='Actual Valuation $v_2(denom)$', markersize=4)
        
        plt.title("Zudilin Row: 2-adic Denominator Growth", fontweight='bold')
        plt.xlabel("Index m")
        plt.ylabel("Power of 2")
        plt.grid(True, linestyle=':', alpha=0.7)
        plt.legend()
        plt.tight_layout()
        plt.savefig("p2_valuation_bound.png", dpi=300, bbox_inches='tight')
        plt.close()
        print("Saved: p2_valuation_bound.png")

    # =========================================================================
    # PLOT 2: Digits Correct vs Index n
    # =========================================================================
    plt.figure(figsize=(7, 5))
    max_n = 25

    if zudilin_data:
        z_n = [n for n in zudilin_data['n'] if n <= max_n]
        z_dig = zudilin_data['digits_correct'][:len(z_n)]
        plt.plot(z_n, z_dig, 'go-', label='Zudilin (at 3n)', markersize=5)
        
    if nesterenko_data:
        plt.plot(nesterenko_data['n'], nesterenko_data['digits_correct'], 
                'bo-', label='Nesterenko (Row 2)', markersize=5)
        
    if combined_data:
        plt.plot(combined_data['n'], combined_data['best_digits_correct'], 
                'mo-', label='Combined Lattice Row', markersize=5)

    plt.title("Approximation Accuracy vs. Index n", fontweight='bold')
    plt.xlabel("Sequence Index n")
    plt.ylabel("Decimal Digits Correct")
    plt.xlim(0, max_n)
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig("accuracy_vs_index.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: accuracy_vs_index.png")

    # =========================================================================
    # PLOT 3: Digits Correct vs Denominator Bit-Length (The LCM Gain)
    # =========================================================================
    plt.figure(figsize=(7, 5))
    
    if zudilin_data:
        plt.plot(zudilin_data['X_n_bit_length'], zudilin_data['digits_correct'], 
                'go-', label='Zudilin Row', markersize=5)
        
    if nesterenko_data:
        plt.plot(nesterenko_data['V_n_bit_length'], nesterenko_data['digits_correct'], 
                'bo-', label='Nesterenko Row', markersize=5)
        
    if combined_data:
        plt.plot(combined_data['best_q_n_bit_length'], combined_data['best_digits_correct'], 
                'mo-', label='Combined Row (LCM-Square)', markersize=5)

    plt.title("Efficiency: Digits Correct per Bit", fontweight='bold')
    plt.xlabel("Denominator Bit-Length")
    plt.ylabel("Decimal Digits Correct")
    plt.grid(True, linestyle=':', alpha=0.7)
    plt.legend()
    plt.tight_layout()
    plt.savefig("efficiency_digits_per_bit.png", dpi=300, bbox_inches='tight')
    plt.close()
    print("Saved: efficiency_digits_per_bit.png")

if __name__ == "__main__":
    main()
