import math

def orthogonal_projection(basis_vec, vec):
    norm_basis = math.sqrt(sum([x**2 for x in basis_vec]))
    proj_coeff = sum([basis_vec[i] * vec[i] for i in range(len(basis_vec))]) / (norm_basis ** 2)
    return [proj_coeff * basis_vec[i] for i in range(len(basis_vec))]

A = [[1, 0.00000001, 0, 0],
     [1, 0, 0.00000001, 0],
     [1, 0, 0, 0.00000001]]

basis = []
orthogonal_basis_classical = []

for j in range(len(A)):
    basis.append(A[j][:])

    for i in range(j):
        proj = orthogonal_projection(basis[i], A[j])
        basis[j] = [basis[j][k] - proj[k] for k in range(len(A[j]))]

for j in range(len(basis)):
    norm_basis_j = math.sqrt(sum([basis[j][i]**2 for i in range(len(basis[j]))]))
    orthogonal_basis_classical.append([basis[j][k] / norm_basis_j for k in range(len(basis[j]))])


def vector_norm(v):
    return math.sqrt(sum([x ** 2 for x in v]))

def dot_product(v1, v2):
    return sum([v1[i] * v2[i] for i in range(len(v1))])

def modified_gram_schmidt(A):
    num_vectors = len(A)
    Q = [A[j][:] for j in range(num_vectors)]

    for j in range(num_vectors):
        norm_q_j = vector_norm(Q[j])
        if norm_q_j == 0:
            raise ValueError("error")

        Q[j] = [Q[j][i] / norm_q_j for i in range(len(Q[j]))]

        for k in range(j + 1, num_vectors):
            projection_coeff = dot_product(Q[k], Q[j])
            Q[k] = [Q[k][i] - projection_coeff * Q[j][i] for i in range(len(Q[k]))]

    return Q

Q = modified_gram_schmidt(A)

def matrix_product(M1, M2):
    result = [[0] * len(M2[0]) for _ in range(len(M1))]
    for i in range(len(M1)):
        for j in range(len(M2[0])):
            result[i][j] = sum(M1[i][k] * M2[k][j] for k in range(len(M2)))
    return result

def transpose(matrix):
    return [[matrix[j][i] for j in range(len(matrix))] for i in range(len(matrix[0]))]


def create_matrix(orthogonal_basis):
    return [orthogonal_basis[i][:] for i in range(len(orthogonal_basis))]

classical_matrix = create_matrix(orthogonal_basis_classical)
modified_matrix = create_matrix(Q)

classical_transpose = transpose(classical_matrix)
modified_transpose = transpose(modified_matrix)

classical_product = matrix_product(classical_matrix, classical_transpose)
modified_product = matrix_product(modified_matrix, modified_transpose)


print("Classical Gram-Schmidt Matrix:")
for row in classical_matrix:
    print(row)

print("\nClassical Gram-Schmidt Transpose:")
for row in classical_transpose:
    print(row)

print("\nClassical Gram-Schmidt Product with its Transpose:")
for row in classical_product:
    print(row)

print("\nModified Gram-Schmidt Matrix:")
for row in modified_matrix:
    print(row)

print("\nModified Gram-Schmidt Transpose:")
for row in modified_transpose:
    print(row)

print("\nModified Gram-Schmidt Product with its Transpose:")
for row in modified_product:
    print(row)

# https://laurenthoeltgen.name/post/gram-schmidt/