from my_functions import *

normal=Prob_a_to_b_General("mu","mu","normal")
anti=Prob_a_to_b_General("mu","mu","anti")

print("Prob of normal:",normal)
print("Prob of anti:",anti)

print("Difference of the two:",simplify( normal-anti))
