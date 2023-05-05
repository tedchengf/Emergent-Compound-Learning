import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
from Models import *

# Global Variables
global elements, products, emergent_compound_products
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

def main():
	# Define Ground Truth Model
	ground_truth = Compound_Model(elements, products, emergent_compound_products)

	ps = Partition_Space(ground_truth, constrain_compounds=constrain_compounds)

	# data = np.loadtxt("simulation_output.txt", dtype = float, skiprows = 1)
	# all_scores = get_problem_score(ground_truth.compounds,ps.hypotheses_incidence)
	
	# fig, axs = plt.subplots(4, 1, sharex = True, figsize = (16, 12))
	# fig.subplots_adjust(hspace=0.1)
	# matplotlib.rcParams.update({'font.size': 2})
	# colors = plt.rcParams["axes.prop_cycle"].by_key()['color'][:4]
	# axs[0].scatter(data[:, 0].flatten(), data[:, 1].flatten(), color = colors[0], alpha = 0.15, s = 30)
	# axs[0].set_ylabel("Trial Num")
	# axs[1].scatter(data[:, 0].flatten(), data[:, 2].flatten(), color = colors[1], alpha = 0.15, s = 30)
	# axs[1].set_ylabel("Expected IG")
	# axs[2].scatter(data[:, 0].flatten(), data[:, 3].flatten(), color = colors[2], alpha = 0.15, s = 30)
	# axs[2].set_ylabel("Actual IG")
	# sns.kdeplot(all_scores, ax = axs[3], fill = True, color = "grey")
	# sns.set_style("darkgrid")
	# axs[3].set_ylabel("Density")
	# axs[3].set_xlabel("Problem Complexity")
	# fig.savefig("Scatters.png", format = "png", dpi = 1000, transparent = True)
	# return

	# for ind in range(4):
	# 	axs[ind].scatter(score, ys[ind], color = colors[ind], alpha = 0.7, s = 12)
	# 	axs[ind].set_ylabel(var_names[ind])
	# axs[-1].set_xlabel("Problem Complexity")
	# fig.savefig("Scatters.png", format = "png", dpi = 1000, transparent = True)	

	# fig = plt.figure(figsize = (16, 3))
	# ax = fig.gca()
	# sns.kdeplot(all_scores, ax = ax, fill = True, color = "grey")
	# sns.set_style("darkgrid")
	# # ax.set_xticklabels([])
	# # ax.set_xticks(np.arange(0.2, 1.1, 0.1))
	# fig.savefig("Score Distribution.png", format = "png", dpi = 1000, transparent = True)
	# return

	# curr_gamma = lambda x, y: gamma_prior(x,y,scaling = 5)
	# sc, tn, eig, aig = simulate_problems(ps, ground_truth.compounds, curr_gamma, simple_llh, entropy, repeat = 5, sample = 200)
	# with open("simulation_output.txt", "w") as outfile:
	# 	outfile.write("Scores\tTrial Num\tEIG\tAIG\n")
	# 	for ind in range(len(sc)):
	# 		outfile.write(str(sc[ind]))
	# 		outfile.write("\t")
	# 		outfile.write(str(tn[ind]))
	# 		outfile.write("\t")
	# 		outfile.write(str(eig[ind]))
	# 		outfile.write("\t")
	# 		outfile.write(str(aig[ind]))
	# 		outfile.write("\t")
	# 		outfile.write("\n")

	# stacked_plot(sc, [tn, eig, aig], ["Trial Num", "EIG","AIG"], (16, 9))
	# simple_plot(sc, tn, (12.4, 4.8), "orange", "Trial Num")
	# simple_plot(sc, eig, (12.4, 4.8), "blue", "Expected Information Gain")
	# simple_plot(sc, eig, (12.4, 4.8), "pink", "Actual Information Gain")
	
	# return

	# Define a prior distribution function
	curr_gamma = lambda x, y: gamma_prior(x,y,scaling = 0.5)
	# Define an active learning agent
	ac_agent = Active_Learning(ground_truth, constrain_compounds, uniform_prior, simple_llh)
	# Perform cctive learning
	best_model, model_prob, choice_history, posterior_entropy, expected_IG, actual_IG = ac_agent.active_learning()
	print(choice_history)
	print(posterior_entropy)
	print(expected_IG)
	print(actual_IG)

def stacked_plot(score, ys, var_names, figsize):
	fig, axs = plt.subplots(len(ys), 1, sharex = True, figsize = figsize)
	fig.subplots_adjust(hspace=0.1)
	matplotlib.rcParams.update({'font.size': 2})
	colors = plt.rcParams["axes.prop_cycle"].by_key()['color'][:len(ys)]
	for ind in range(len(ys)):
		axs[ind].scatter(score, ys[ind], color = colors[ind], alpha = 0.7, s = 12)
		axs[ind].set_ylabel(var_names[ind])
	axs[-1].set_xlabel("Problem Complexity")
	fig.savefig("Scatters.png", format = "png", dpi = 1000, transparent = True)

def simple_plot(x, y, figsize, color, var_name):
	fig = plt.figure(figsize = figsize)
	ax = fig.gca()
	ax.scatter(x, y, color = color)
	ax.set_xlabel("score")
	ax.set_ylabel(var_name)
	fig.savefig(var_name + ".png", format = "png", dpi = 1000, transparent = True)

def simulate_problems(partition_space, compounds, prior_func, llh_func, ig_func, constrain = 5, repeat = 5, sample = 200):
	simulate_prods = 1 + np.arange(partition_space.hypotheses_incidence.shape[1])
	start = np.vstack(([partition_space.hypotheses_incidence[0]], [partition_space.hypotheses_incidence[-1]]))
	rest_idx = np.random.choice(np.arange(partition_space.hypotheses_incidence.shape[0] - 2, dtype = int), sample - 2, replace = False)
	rest = partition_space.hypotheses_incidence[1:-1][rest_idx]
	all_tests = np.vstack((start, rest))
	scores = get_problem_score(compounds, all_tests)
	all_scores = []
	trial_num = []
	averaged_EIG = []
	averaged_AIG = []

	curr_ind = 0
	for test in tqdm(all_tests):
		ground_truth_dict = {}
		constrain_dict = {}
		curr_partition = compounds[test.astype(bool)]
		curr_prods = list(simulate_prods[test.astype(bool)])
		for ind in range(len(curr_prods)):
			comp = curr_partition[ind]
			prod = curr_prods[ind]
			ground_truth_dict.update({comp: prod})
			if ind < constrain:
				constrain_dict.update({comp: prod})

		curr_ground_truth = Compound_Model(elements, curr_prods, ground_truth_dict)
		model_trial_num = []
		model_averaged_EIG = []
		model_averaged_AIG = []
		for ind in range(repeat):
			curr_ac = Active_Learning(curr_ground_truth, constrain_dict, prior_func, llh_func, ig_func)
			best_model, model_prob, choice_history, posterior_entropy, expected_IG, actual_IG = curr_ac.active_learning(return_n = 1, verbose = False)
			all_scores.append(scores[curr_ind])
			trial_num.append(len(choice_history))
			averaged_EIG.append(np.mean(expected_IG[1:]))
			averaged_AIG.append(np.mean(actual_IG[1:]))
		
		# 	model_trial_num.append(len(choice_history))
		# 	model_averaged_EIG.append(np.mean(expected_IG[1:]))
		# 	model_averaged_AIG.append(np.mean(actual_IG[1:]))
		# trial_num.append(np.mean(model_trial_num))
		# averaged_EIG.append(np.mean(model_averaged_EIG))
		# averaged_AIG.append(np.mean(model_averaged_AIG))
		
		curr_ind += 1
	
	return all_scores, trial_num, averaged_EIG, averaged_AIG

def get_problem_score(compounds, incidences, alpha = 1, beta = 1):
	all_card = len(compounds)
	all_size = 0
	for c in compounds: all_size += len(c)
	all_x = np.empty(incidences.shape[0])
	for ind, inc in enumerate(incidences):
		curr_card = np.sum(inc)
		curr_compounds = compounds[inc.astype(bool)]
		curr_size = 0
		for c in curr_compounds: curr_size += len(c)
		# Here x is in the range of [0,1]
		all_x[ind] = (alpha*(curr_card/all_card) + beta*(curr_size/all_size))/2
	return all_x

if __name__ == "__main__":
	main()