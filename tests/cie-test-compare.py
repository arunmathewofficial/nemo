import numpy as np
import sys


def compare_vectors(a, b, rel_tol=3e-5, abs_tol=3e-11):
    """
  Compare two vectors with relative and absolute tolerance
  Returns True if all elements are within tolerance, False otherwise
  """
    if len(a) != len(b):
        raise ValueError("Vectors must have the same length")

    # Compare each element
    for i in range(len(a)):
        diff = abs(a[i] - b[i])
        abs_sum = abs(a[i]) + abs(b[i])
        rel_diff = diff / abs_sum if abs_sum != 0 else 0
        if diff > abs_tol and rel_diff > rel_tol:
            print("i=", i, "a=", a[i], "b=", b[i])
            return False
    return True


def main():
    # Check command line arguments
    if len(sys.argv) < 3:
        print(
            "Usage: python compare.py <output_file> <comparison_file> [<rel-tol>] [<abs-tol>]"
        )
        return 1

    output_file = sys.argv[1]
    comparison_file = sys.argv[2]
    if (len(sys.argv) >= 4):
        reltol = sys.argv[3]
        reltol = float(reltol)
    else:
        reltol = 3.0e-5

    if (len(sys.argv) == 5):
        abstol = sys.argv[4]
        abstol = float(abstol)
    else:
        abstol = 1.0e-14

    # Load the data files
    try:
        output_data = np.loadtxt(output_file, skiprows=1)
        comparison_data = np.loadtxt(comparison_file, skiprows=1)

        # Check if files have the same number of rows
        if output_data.shape[0] != comparison_data.shape[0]:
            print(
                f"Warning: Different number of rows - output: {output_data.shape[0]}, comparison: {comparison_data.shape[0]}"
            )
            min_rows = min(output_data.shape[0], comparison_data.shape[0])
            return 1
        else:
            min_rows = output_data.shape[0]

        # Compare each row
        all_match = True
        for i in range(min_rows):
            # Compare all data on a row
            output_row = output_data[i]
            comparison_row = comparison_data[i]

            if not compare_vectors(
                    output_row, comparison_row, rel_tol=reltol,
                    abs_tol=abstol):
                print(f"Row {i} differs")
                print("op", output_row)
                print("cm", comparison_row)
                all_match = False
            else:
                print(f"Row {i} matches")

        if all_match:
            print("All rows match within tolerance")
            return 0
        else:
            print("Some rows differ")
            return 1

    except FileNotFoundError as e:
        print(f"File not found: {e}")
        return 1
    except Exception as e:
        print(f"Error during comparison: {e}")
        return 1


if __name__ == "__main__":
    exit(main())
