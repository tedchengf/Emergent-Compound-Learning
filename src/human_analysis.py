import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from Models import *
from scipy.special import softmax

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

global subnames, sub_pass
sub_names = ["fc01_processed.csv",
	     	 "fc02_processed.csv",
			 "fc03_processed.csv",
			 "fc04_processed.csv",
			 "fc05_processed.csv",
			 "hl02_processed.csv",
			 "hl03_processed.csv",
			 "sl01_processed.csv",
			 "sl02_processed.csv",
			 "sl03_processed.csv",
			 "sl05_processed.csv",
			 "sy01_processed.csv",
			 "sy02_processed.csv",
			 "sy03_processed.csv",]
sub_pass = np.array([1,1,1,1,1,0,1,0,1,0,1,1,0,0], dtype = bool)

def main():
	# Define Ground Truth Model
	ground_truth = Compound_Model(elements, products, emergent_compound_products)
	# curr_gamma = lambda x, y: gamma_prior(x,y,scaling = 0.5)
	# # Define an active learning agent
	# ac_agent = Active_Learning(ground_truth, constrain_compounds,
	# curr_gamma, simple_llh)
	# sub_compounds, rt1, rt2 = read_choice_data("./Human_Data/data_game.csv")

	for sub_ind in tqdm(range(len(sub_names))):
		# if sub_ind == 2: continue
		# if sub_ind == 10: continue
		# if sub_ind == 13: continue
		sub_name = sub_names[sub_ind]
		analyze_subject(ground_truth, sub_name, "./Subject_Analysis/")

	# # Main Loop
	# best_IG, sub_IG, random_IG, worst_IG, choice_entrophy, suprisal, compound_prob, random_prob, selection_prob = [],[],[],[],[],[],[],[],[]
	# for ind in range(len(sub_compounds)):
	# 	best_compound, IG, compound_ranks, compounds, dist = ac_agent.choose_compound()
	# 	choice_entrophy.append(entropy(normalize_arr(dist)))
	# 	BIG, SIG, RIG, WIG,s_prob, r_prob = choice_entropy_deviation(best_compound, sub_compounds[ind], compounds, dist)
	# 	compound_prob.append(ac_agent.partition_space.compound_prob(sub_compounds[ind]))
	# 	best_IG.append(BIG)
	# 	sub_IG.append(SIG)
	# 	random_IG.append(RIG)
	# 	worst_IG.append(WIG)
	# 	selection_prob.append(s_prob)
	# 	random_prob.append(r_prob)
	# 	# Always remember to advance the model by putting down the actual
	# 	# compound the subject chooses
	# 	models, posterior, AIG = ac_agent.learning_episode(sub_compounds[ind])
	# 	suprisal.append(abs(AIG - SIG))
	
	# # Analysis 1: Sequential Information Gain
	# fig = plt.figure(figsize=(16,4))
	# ax = fig.gca()
	# ax.plot(best_IG, label = "Best Choice IG")
	# ax.plot(sub_IG, label = "Subject Choice IG")
	# ax.plot(random_IG, label = "Random Choice IG")
	# ax.plot(worst_IG, label = "Worst Choice IG")
	# ax.set_xlabel("Trials")
	# ax.set_ylabel("Gain in Shannon Entropy")
	# fig.legend()
	# fig.savefig("Sequential IG.png", format = "png", dpi = 1000, transparent=True)

	# # Analysis 1.5: probability
	# fig = plt.figure(figsize=(16,4))
	# ax = fig.gca()
	# ax.plot(selection_prob, label = "Subject Choice probability")
	# ax.plot(random_prob, label = "Random Choice probability")
	# ax.set_xlabel("Trials")
	# ax.set_ylabel("Model Probability")
	# fig.legend()
	# fig.savefig("Sequential Prob.png", format = "png", dpi = 1000, transparent=True)

	# # Analysis 2: Linear Testing Strategy
	# fig = plt.figure(figsize=(16,4))
	# ax = fig.gca()
	# ax.plot(compound_prob, label = "Probability of Compound being and Emergent Compound")
	# ax.set_xlabel("Trials")
	# ax.set_ylabel("Probability")
	# fig.legend()
	# fig.savefig("Linear Testing.png", format = "png", dpi = 1000, transparent=True)

	# # Analysis 2: Reaction time and uncerstainty
	# fig, axs = plt.subplots(2, 1, sharex = True, figsize = (16, 6))
	# fig.subplots_adjust(hspace=0.1)
	# axs[0].plot(normalize_arr(rt1), label = "Choice Reaction Time")
	# axs[0].plot(normalize_arr(choice_entrophy), label = "Choice Entrophy")
	# axs[1].plot(normalize_arr(rt2), label = "After Choice Reaction Time")
	# axs[1].plot(normalize_arr(suprisal), label = "Surprisal")
	# axs[0].set_ylabel("Normalized Reaction Time / Entrophy")
	# axs[1].set_ylabel("Normalized Reaction Time / Entrophy")
	# axs[1].set_xlabel("Trials")
	# fig.legend()
	# fig.savefig("Reaction Times.png", format = "png", dpi = 1000, transparent=True)
	return

##### TODO: why the uniform dist does not change the results?
def analyze_subject(ground_truth, subfile, prefix = ""):
	curr_gamma = lambda x, y: gamma_prior(x,y,scaling = 0.5)
	# Define an active learning agent
	ac_agent = Active_Learning(ground_truth, constrain_compounds,
	curr_gamma, simple_llh)
	sub_compounds, elements, rt1, rt2 = read_choice_data("./Subject_Data/"+subfile)

	# Main Loop
	best_IG, sub_IG, random_IG, worst_IG, choice_entrophy, suprisal, compound_prob, best_prob, random_prob, selection_prob, worst_prob, actual_elements = [],[],[],[],[],[],[],[],[],[],[],[]
	for ind in range(len(sub_compounds)):
		best_compound, IG, compound_ranks, compounds, dist = ac_agent.choose_compound()
		choice_entrophy.append(entropy(normalize_arr(dist)))
		try:
			BIG, SIG, RIG, WIG,b_prob,s_prob,r_prob,w_prob = choice_entropy_deviation(best_compound, sub_compounds[ind], compounds, dist)
		except UnboundLocalError: 
			continue
		actual_elements.append(elements[ind])
		compound_prob.append(ac_agent.partition_space.compound_prob(sub_compounds[ind]))
		best_IG.append(BIG)
		sub_IG.append(SIG)
		random_IG.append(RIG)
		worst_IG.append(WIG)
		best_prob.append(b_prob)
		selection_prob.append(s_prob)
		random_prob.append(r_prob)
		worst_prob.append(w_prob)
		# Always remember to advance the model by putting down the actual
		# compound the subject chooses
		models, posterior, AIG = ac_agent.learning_episode(sub_compounds[ind])
		suprisal.append(abs(AIG - SIG))
	
	matplotlib.rcParams.update({'font.size': 12})

	fig, axs = plt.subplots(3, 1, sharex = False, figsize = (11, 10))
	fig.subplots_adjust(hspace=0.1)

	# Analysis 1: Sequential Information Gain
	axs[0].plot(best_IG, label = "Optimal Choice Utility")
	axs[0].plot(sub_IG, label = "Subject Choice Utility")
	axs[0].plot(random_IG, label = "Random Choice Utility")
	axs[0].plot(worst_IG, label = "Worst Choice Utility")
	axs[0].set_ylabel("Change in Shannon Entropy")
	axs[0].xaxis.set_label_position('top')
	axs[0].tick_params(top=True, labeltop=True, bottom=False, labelbottom=False)
	axs[0].set_xticks(list(np.arange(len(best_IG), dtype = int)))
	axs[0].set_xticklabels(np.arange(len(best_IG), dtype = int).astype(str))
	axs[0].set_xlabel("Trials / Subject's choice of the compound to test")
	axs[0].legend()

	# Analysis 1.5: probability
	axs[1].plot(best_prob,label = "Optimal Choice probability")
	axs[1].plot(selection_prob, label = "Subject Choice probability")
	axs[1].plot(random_prob, label = "Random Choice probability")
	axs[1].plot(worst_prob, label = "Worst Choice probability")
	axs[1].set_ylabel("Model Probability")
	axs[1].set_xticks(list(np.arange(len(best_IG), dtype = int)))
	axs[1].set_xticklabels([])
	axs[1].legend()

	# Analysis 2: Linear Testing Strategy
	axs[2].plot(compound_prob, label = "Probability of Compound being an Emergent Compound", color = "orange")
	axs[2].hlines(0, xmin = 0, xmax = len(best_IG) - 1, colors = "red", label = "0 Probability (Violation of Linearity)")
	axs[2].hlines(1, xmin = 0, xmax = len(best_IG) - 1, colors = "red")
	# axs[2].set_xlabel("Trials (Subject's choice of the compound to test)")
	axs[2].set_xticks(np.arange(len(best_IG), dtype = int))
	axs[2].set_xticklabels(actual_elements, rotation=30)
	axs[2].set_ylabel("Model Probability")
	axs[2].legend()
	fig.savefig(prefix + subfile[:4] + "_Results.png", format = "png", dpi = 500, transparent=True)
	plt.close()

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
	prob = softmax(dist)
	# print(dist)
	# print(prob)
	best_IG = np.argwhere(compounds == best_compound)
	for ind, c in enumerate(compounds):
		c = tuple(c)
		if c == best_compound: best_ind = ind
		if c == subject_compound: sub_ind = ind
	best_IG = dist[best_ind]
	try:
		sub_IG = dist[sub_ind]
	except UnboundLocalError:
		print(subject_compound)
	random_IG = np.mean(dist)
	worst_IG = np.amin(dist)
	best_prob = prob[best_ind]
	selection_prob = prob[sub_ind]
	random_prob = 1/len(compounds)
	worst_prob = np.amin(prob)
	return best_IG, sub_IG, random_IG, worst_IG, best_prob, selection_prob, random_prob, worst_prob

def read_choice_data(fname):
	data = np.loadtxt(fname, skiprows=1, dtype=str, delimiter = ",")
	# print(data)
	# print(data.shape)
	choices = data[:,0]
	elements = []
	compounds = []
	for l in choices:
		elems = l.strip("\"").split()
		elements.append(" ".join(elems))
		curr_comp = [inputs_to_elements[e] for e in elems]
		compounds.append(tuple(sorted(curr_comp)))
	rt1 = data[:,1].flatten().astype(int)
	rt2 = data[:,2].flatten().astype(int)
	return compounds, elements, rt1, rt2

if __name__ == "__main__":
	main()