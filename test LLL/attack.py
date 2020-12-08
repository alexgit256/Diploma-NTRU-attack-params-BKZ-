#!/usr/bin/env g6k

from fpylll import BKZ as fplll_bkz, GSO, IntegerMatrix, LLL
import time

def read_ntru_from_file(filename):

    data = open(filename, "r").readlines()
    q = int(data[0])
    H = eval(",".join([s_.replace('\n','').replace(" ", ", ") for s_ in data[1 :]]))
    H = IntegerMatrix.from_matrix(H)
    return H, q

def ntru_plain_hybrid_basis(A, k, q, time_flag=True):
	"""
		Construct ntru lattice basis
	"""
	
	if time_flag:
		T=time.perf_counter()
	#print(type(A))
	n = A.ncols
	print("k=", k)
	
	B = IntegerMatrix(2*k,2*k)
	
	#generate Identity matrix of side k x k in upper left corner
	for i in range(0,k):
		B[i,i]=1
	#generate q*Identity matrix of side k x k in lower right corner
	for i in range(k,2*k):
		B[i,i]=q

	#generate upper right block of size kxk
	for j in range(0,k):
		for l in range(0,k):
			B[j,k+l]=A[n-k+j,n-k+l]
	#print(B)
	B = LLL.reduction(B)
	if time_flag:
		print(time.perf_counter()-T, "sec for LLL")
	
	Bg=B
	return B, Bg
	
def write_basis(B, q, filename):
	f = open(filename, 'w')
	f.write(str(q)+'\n')
	for t in B:
		f.write( str(t).replace(',','') +'\n')
	f.close()

A, q= read_ntru_from_file("///home/alex/proj/g6k/g6k/diploma/ntru_n_64.txt")
k=62
print(q)
#print(A)
B, Bg = ntru_plain_hybrid_basis(A, k, q)
filename='ntru_k_'+str(k)+'_LLL.txt'
write_basis(B,q,filename)
