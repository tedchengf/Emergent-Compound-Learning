import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from Models import *

# Global Variables
global elements, products, emergent_compound_products, inputs_to_elements
elements = ["A","B","C","D"]
products = ["0", "♁", "♆", "🜺", "🜍", "♄", "☽", "☉"]
emergent_compound_products = {
	(): products[0],
	("A",): products[1],
	("B",): products[2],
	("C",): products[3],
	("D",): products[4],
	("A","B"): products[5],
	("A", "B", "D"): products[6],
	("B", "C", "D"): products[7]
}
constrain_compounds = {
	(): "0",
	("A",): emergent_compound_products[("A",)],
	("B",): emergent_compound_products[("B",)],
	("C",): emergent_compound_products[("C",)],
	("D",): emergent_compound_products[("D",)]
}
inputs_to_elements = {
	"Air": "A",
	"Earth": "B",
	"Fire": "C",
	"Water": "D"
}

def main():
	# Define Ground Truth Model
	ground_truth = Compound_Model(elements, products, emergent_compound_products)
	curr_gamma = lambda x, y: gamma_prior(x,y,scaling = 0.5)
	# Define an active learning agent
	ac_agent = Active_Learning(ground_truth, constrain_compounds,
	curr_gamma, simple_llh)
	sub_compounds, rt1, rt2 = read_choice_data("./Human_Data/data_game.csv")


	# Main Loop
	best_IG, sub_IG, random_IG, worst_IG, choice_entrophy, suprisal = [],[],[],[], [], []
	for ind in range(len(sub_compounds)):
		best_compound, IG, compound_ranks, compounds, dist = ac_agent.choose_compound()
		choice_entrophy.append(entropy(normalize_arr(dist)))
		BIG, SIG, RIG, WIG = choice_entropy_deviation(best_compound, sub_compounds[ind], compounds, dist)
		best_IG.append(BIG)
		sub_IG.append(SIG)
		random_IG.append(RIG)
		worst_IG.append(WIG)
		# Always remember to advance the model by putting down the actual
		# compound the subject chooses
		models, posterior, AIG = ac_agent.learning_episode(sub_compounds[ind])
		suprisal.append(abs(AIG - SIG))
	
	# Analysis 1: Sequential Information Gain
	fig = plt.figure(figsize=(16,3))
	ax = fig.gca()
	ax.plot(best_IG, label = "Best Choice IG")
	ax.plot(sub_IG, label = "Subject Choice IG")
	ax.plot(random_IG, label = "Random Choice IG")
	ax.plot(worst_IG, label = "Worst Choice IG")
	ax.set_xlabel("Trials")
	ax.set_ylabel("Gain in Shannon Entropy")
	fig.legend()
	fig.savefig("Sequential IG.png", format = "png", dpi = 1000, transparent=True)
	
	# Analysis 2: Reaction time and uncerstainty
	fig, axs = plt.subplots(2, 1, sharex = True, figsize = (16, 6))
	fig.subplots_adjust(hspace=0.1)
	axs[0].plot(normalize_arr(rt1), label = "Choice Reaction Time")
	axs[0].plot(normalize_arr(choice_entrophy), label = "Choice Entrophy")
	axs[1].plot(normalize_arr(rt2), label = "After Choice Reaction Time")
	axs[1].plot(normalize_arr(suprisal), label = "Surprisal")
	axs[0].set_ylabel("Normalized Reaction Time / Entrophy")
	axs[1].set_ylabel("Normalized Reaction Time / Entrophy")
	axs[1].set_xlabel("Trials")
	fig.legend()
	fig.savefig("Reaction Times.png", format = "png", dpi = 1000, transparent=True)
	return

def normalize_arr(x):
	x = np.array(x)
	return x / x.sum()

def simple_plot(x,figsize, color, var_name):
	fig = plt.figure(figsize = figsize)
	ax = fig.gca()
	ax.scatter(x, y, color = color)
	ax.set_xlabel("score")
	ax.set_ylabel(var_name)
	fig.savefig(var_name + ".png", format = "png", dpi = 1000, transparent = True)

def choice_entropy_deviation(best_compound, subject_compound, compounds, dist):
	best_IG = np.argwhere(compounds == best_compound)
	for ind, c in enumerate(compounds):
		c = tuple(c)
		if c == best_compound: best_ind = ind
		if c == subject_compound: sub_ind = ind
	best_IG = dist[best_ind]
	sub_IG = dist[sub_ind]
	random_IG = np.mean(dist)
	worst_IG = np.amin(dist)
	return best_IG, sub_IG, random_IG, worst_IG

def read_choice_data(fname):
	data = np.loadtxt(fname, skiprows=1, dtype=str, delimiter = ",")
	print(data)
	print(data.shape)
	choices = data[:,0]
	compounds = []
	for l in choices:
		elems = l.strip("\"").split()
		curr_comp = [inputs_to_elements[e] for e in elems]
		compounds.append(tuple(sorted(curr_comp)))
	rt1 = data[:,1].flatten().astype(int)
	rt2 = data[:,2].flatten().astype(int)
	return compounds, rt1, rt2

if __name__ == "__main__":
	main()