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
        abs_sum = 0.5 * (abs(a[i]) + abs(b[i]))
        rel_diff = diff / abs_sum if abs_sum != 0 else 0
        if diff > abs_tol and rel_diff > rel_tol:
            print("i=", i, "a=", a[i], "b=", b[i], "reldiff:", rel_diff,
                  "absdiff:", diff)
            return False
    return True


def main():
    # Check command line arguments
    if len(sys.argv) < 3:
        print(
            "Usage: python compare.py <newdata_file> <refdata_file> [<rel-tol>] [<abs-tol>]"
        )
        return 1

    newdata_file = sys.argv[1]
    refdata_file = sys.argv[2]
    if (len(sys.argv) >= 4):
        reltol = sys.argv[3]
        reltol = float(reltol)
    else:
        reltol = 3.0e-5

    if (len(sys.argv) == 5):
        abstol = sys.argv[4]
        abstol = float(abstol)
    else:
        abstol = 1.0e-10

    # Load the data files
    try:
        print("Reading new data:", newdata_file)
        print("Reading ref data:", refdata_file)
        newdata_data = np.loadtxt(newdata_file, skiprows=1)
        refdata_data = np.loadtxt(refdata_file, skiprows=1)

        # Check if files have the same number of rows
        if newdata_data.shape[0] != refdata_data.shape[0]:
            print(
                f"Warning: Different number of rows - output: {newdata_data.shape[0]}, comparison: {refdata_data.shape[0]}"
            )
            min_rows = min(newdata_data.shape[0], refdata_data.shape[0])
            return 1
        else:
            min_rows = newdata_data.shape[0]

        # Compare each row
        all_match = True
        for i in range(min_rows):
            # Compare all data on a row
            newdata_row = newdata_data[i]
            refdata_row = refdata_data[i]

            if not compare_vectors(
                    newdata_row, refdata_row, rel_tol=reltol, abs_tol=abstol):
                print(f"Row {i} differs")
                print("op", newdata_row)
                print("cm", refdata_row)
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
