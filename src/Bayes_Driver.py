import numpy as np
from Models import *

# Global Variables
global elements, emergent_compound_rewards, constrain_compounds
elements = ["A","B","C","D"]
emergent_compound_rewards = {
	(): 0,
	("A",): np.random.normal(),
	("B",): np.random.normal(),
	("C",): np.random.normal(),
	("D",): np.random.normal(),
	("A","B"): 7,
	("B", "D"): 12,
	("A", "C", "D"): -3
}
constrain_compounds = {
	(): 0,
	("A",): emergent_compound_rewards[("A",)],
	("B",): emergent_compound_rewards[("B",)],
	("C",): emergent_compound_rewards[("C",)],
	("D",): emergent_compound_rewards[("D",)]
}
# Random Seed
np.random.seed(2023)

def main():
	ground_truth = Compound_Model(elements, emergent_compound_rewards)
	ps = Partition_Space(ground_truth, constrain_compounds)
	# print(simple_bayes.hypotheses.shape)
	# print(simple_bayes.hypotheses)
	return

if __name__ == "__main__":
	main()