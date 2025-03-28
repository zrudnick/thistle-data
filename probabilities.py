
# input data and count character at each position
import sys

def find_most_likely_chars(ifn):

    with open(ifn, "r") as data:

        # A, C, G, T
        mat = [[0, 0, 0, 0] for i in range(203)]
        for line in data:
            if line[0].isalpha():
                for j in range(len(line[:-1])):
                    if line[j] == 'A':
                        mat[j][0] += 1
                    elif line[j] == 'C':
                        mat[j][1] += 1
                    elif line[j] == 'G':
                        mat[j][2] += 1
                    elif line[j] == 'T':
                        mat[j][3] += 1

        for i in range(len(mat)):
            total = 0
            for j in range(4): total += mat[i][j] 
            for j in range(4): mat[i][j] = mat[i][j] / total 

        chars = []
        for i in range(len(mat)):
            max_p = max(mat[i][0], mat[i][1], mat[i][2], mat[i][3])
            j = mat[i].index(max_p)
            if j == 0:
                chars.append('A')
            elif j == 1:
                chars.append('C')
            elif j == 2:
                chars.append('G')
            elif j == 3:
                chars.append('T')

    print(chars[90:110])
    return chars

def main(ifn, ofn):
    chars = find_most_likely_chars(ifn)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

# My Sequence:    C(GCC)GCCGCCATGG
# Kozak Sequence:  (gcc)gccRccAUGG
