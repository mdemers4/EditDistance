import sys
import pprint
from Bio import pairwise2 
from Bio.pairwise2 import format_alignment

def open_txt(argv, n):
	"""
	Opens the text files of interest from arguments in the terminal
	"""
	txt = open(sys.argv[n],"r").read()
	return txt 

def edit_dist(txt1,txt2):
	"""
	Creates a grid of i and j based on the lengths of the txt which
	will calculate the edit distance  
	"""
	first_len = len(txt1)
	second_len = len(txt2)

	DP = [[0]*(second_len+1) for _ in range(first_len+1)]

	for i in range(0,first_len+1): DP[i][0] = i
	#print DP
	for j in range(0,second_len+1): DP[0][j] = j
	#print DP

	for i in range (0,first_len+1):
		for j in range (0,second_len+1):
			if txt1[i-1] == txt2[j-1]:
				DP[i][j] = DP[i-1][j-1]
				

			else:
				DP[i][j] = min(DP[i][j-1],DP[i-1][j],DP[i-1][j-1]) + 1
	return DP[i][j]


def align_sequence(txt1,txt2):
	"""
	accepts two txt files and aligns and prints the 
	sequence based on the its best possible pairwise score.
	"""
	from Bio import pairwise2

	from Bio.SubsMat import MatrixInfo as matlist

	matrix = matlist.blosum62

	txt_list = []

	alns = pairwise2.align.globalms(txt1, txt2, 2, -1, -.1, -1)

	top_aln = alns[0]

	new_txt1, new_txt2, score, begin, end = top_aln


	for i in range(0,len(new_txt1)):
		if new_txt1[i] == new_txt2[i]:
			txt_list.append("|")
		else:
			txt_list.append(" ")
	


	return [new_txt1 , ("".join(map(str, txt_list))) , new_txt2 ]
	



	

if __name__ == "__main__":
	first_text = open_txt(sys.argv, 1) 
	second_text = open_txt(sys.argv, 2)

	print "The edit distance is: %i" % edit_dist(first_text, second_text)

	a = align_sequence(first_text,second_text)


	start=0
	end=59

	while end < (len(a[1])+60):
		for index,i in enumerate(a):
			print a[index][start:end]
		start = end
		end = end + 60



































