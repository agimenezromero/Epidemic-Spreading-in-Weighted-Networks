import numpy as np
import matplotlib.pyplot as plt

from matplotlib import rc

plt.rcParams.update({'font.size': 16})

rc('text', usetex=True)

def theo_lambda_WS(alpha):

	return 8 / (66 * alpha)

def powerlaw():

	plt.figure(figsize=(8, 6))

	data = np.loadtxt("Results/powerlaw.txt")

	t = np.linspace(1, 999, 999)

	plt.plot(t, data[:, 0], ls='-', marker='.', label=r'$\gamma=2$')
	plt.plot(t, data[:, 1], ls='', marker='.', label=r'$\gamma=3$')
	plt.plot(t, data[:, 2], ls='', marker='.', label=r'$\gamma=4$')

	plt.xscale('log')
	plt.yscale('log')

	plt.xlabel(r'$\omega$')
	plt.ylabel(r'$p(\omega)$')

	plt.show()


def alpha_study_SIS_WS(names):

	plt.figure(figsize=(8, 6))

	markers = ['o', '^', 's']

	i = -1

	for name in names:

		i += 1

		data = np.loadtxt("Results/alpha_study_WS_SIS_%s.txt" % name)

		alphas = data[:, 0]
		rho = data[:, 1]

		plt.plot(alphas, rho, ls='-', marker=markers[i], markerfacecolor='None', markersize=10, label=name)

	plt.legend()
	plt.show()

def alpha_study_SIR_WS(names):

	plt.figure(figsize=(8, 6))

	markers = ['o', '^', 's']

	i = -1

	for name in names:

		i += 1

		data = np.loadtxt("Results/alpha_study_WS_SIR_%s.txt" % name)

		alphas = data[:, 0]
		rho = data[:, 1]

		plt.plot(alphas, rho, ls='-', marker=markers[i], markerfacecolor='None', markersize=10, label=name)

	plt.legend()
	plt.show()

def alpha_study_SIR_BA(names):

	plt.figure(figsize=(8, 6))

	markers = ['o', '^', 's']

	i = -1

	for name in names:

		i += 1

		data = np.loadtxt("Results/alpha_study_BA_SIR_%s.txt" % name)

		alphas = data[:, 0]
		rho = data[:, 1]

		plt.plot(alphas, rho, ls='-', marker=markers[i], markerfacecolor='None', markersize=10, label=name)

	plt.legend()
	plt.show()	


def lambda_study_SIS_WS(names):

	plt.figure(figsize=(8, 6))

	markers = ['o', '^', 's']

	i = -1

	for name in names:

		i += 1

		data = np.loadtxt("Results/lambda_study_WS_SIS_%s.txt" % name)

		alphas = data[:, 0]
		rho = data[:, 1]

		plt.plot(alphas, rho, ls='', marker=markers[i], markerfacecolor='None', markersize=10, label=name)

	plt.plot(alphas, theo_lambda_WS(alphas), ls='-', color='k', label='Theory')

	plt.xlim(1, 10)
	plt.ylim(0, 0.13)

	plt.legend()
	plt.show()	

def lambda_study_SIR_WS(names):

	plt.figure(figsize=(8, 6))

	markers = ['o', '^', 's']

	i = -1

	for name in names:

		i += 1

		data = np.loadtxt("Results/lambda_study_WS_SIR_%s.txt" % name)

		alphas = data[:, 0]
		rho = data[:, 1]

		plt.plot(alphas, rho, ls='', marker=markers[i], markerfacecolor='None', markersize=10, label=name)

	plt.plot(alphas, theo_lambda_WS(alphas), ls='-', color='k', label='Theory')

	plt.xlim(1, 10)
	plt.ylim(0, 0.14)

	plt.legend()
	plt.show()	


def dR_lambda(lambdas):

	markers = ['o', '^', 's']

	i = -1

	plt.figure(figsize=(8, 6))

	for lambda_ in lambdas:

		i += 1

		data = np.loadtxt("dR_SIR_lambda_%s.txt" % str(lambda_))

		t = np.linspace(0, len(data), len(data))

		if i == 1:

			plt.plot(t, data, marker=markers[i], label=r'$\lambda_c=%.3f$' % lambda_)

		else:

			plt.plot(t, data, marker=markers[i], label=r'$\lambda=%.3f$' % lambda_)

	plt.plot(np.linspace(0, 5, 100), np.linspace(0.01, 0.01, 100), color='k', linewidth=2)

	plt.xlim(0, 20)

	plt.legend()
	plt.show()


################################################################
#				PAPER 2
#################################################################

def prevalence_t_uniform(ns):

	plt.figure(figsize=(8, 6))

	colors = ['g', 'r', 'b']
	markers = ['s', 'o', '^']

	i = -1

	for n in ns:

		i += 1

		data = np.loadtxt("Results/prevalence_RCM_uniform_n_%s.txt" % str(n))

		t = np.linspace(0, len(data), len(data))

		#plt.plot(t, data, ls='-', color=colors[i], label=r'$n=%i$' % n)
		plt.plot(t[0:-1:20], data[0:-1:20], ls='-', color=colors[i], marker=markers[i], markersize=10, label=r'$n=%i$' % n)

	data = np.loadtxt("Results/random_selection.txt")

	t = np.linspace(0, len(data), len(data))

	plt.plot(t[0:-1:20], data[0:-1:20], ls='--', color='k', lw=2, label='Random selection')

	plt.xlabel(r"$t$")
	plt.ylabel(r"$I(t)$")

	plt.legend(loc='lower right')
	plt.show()

def prevalence_t_powerlaw(gammas):

	plt.figure(figsize=(8, 6))

	colors = ['g', 'r', 'b']
	markers = ['s', 'o', '^']

	i = -1

	for gamma in gammas:

		i += 1

		data = np.loadtxt("Results/prevalence_RCM_powerlaw_gamma_%s.txt" % str(gamma))

		t = np.linspace(0, len(data), len(data))

		#plt.plot(t, data, ls='-', color=colors[i], label=r'$n=%i$' % n)
		plt.plot(t[0:-1:20], data[0:-1:20], ls='-', color=colors[i], marker=markers[i], markersize=10, label=r'$\gamma=%i$' % gamma)


	data = np.loadtxt("Results/random_selection.txt")

	t = np.linspace(0, len(data), len(data))

	plt.plot(t[0:-1:20], data[0:-1:20], ls='--', color='k', lw=2, label='Random selection')

	plt.xlabel(r"$t$")
	plt.ylabel(r"$I(t)$")

	plt.legend(loc='lower right')
	plt.show()

def prevalence_t_uniform_WS(ns):

	plt.figure(figsize=(8, 6))

	colors = ['g', 'r', 'b']
	markers = ['s', 'o', '^']

	i = -1

	for n in ns:

		i += 1

		data = np.loadtxt("Results/prevalence_WS_uniform_n_%s.txt" % str(n))

		t = np.linspace(0, len(data), len(data))

		#plt.plot(t, data, ls='-', color=colors[i], label=r'$n=%i$' % n)
		plt.plot(t[0:-1:20], data[0:-1:20], ls='-', color=colors[i], marker=markers[i], markersize=10, label=r'$n=%i$' % n)

	data = np.loadtxt("Results/prevalence_WS_random.txt")

	t = np.linspace(0, len(data), len(data))

	plt.plot(t[0:-1:20], data[0:-1:20], ls='--', color='k', lw=2, label='Random selection')

	plt.xlabel(r"$t$")
	plt.ylabel(r"$I(t)$")

	plt.legend(loc='lower right')
	plt.show()

def prevalence_t_small_world(ns):

	plt.figure(figsize=(8, 6))

	colors = ['g', 'r', 'b']
	markers = ['s', 'o', '^']

	i = -1

	for n in ns:

		i += 1

		data = np.loadtxt("Results/small_world_lambda_03_n_%s.txt" % str(n))

		t = np.linspace(0, len(data), len(data))

		#plt.plot(t, data, ls='-', color=colors[i], label=r'$n=%i$' % n)
		plt.plot(t[0:-1:100], data[0:-1:100], ls='-', color=colors[i], marker=markers[i], markersize=10, label=r'$n=%i$' % n)

	data = np.loadtxt("Results/small_world_lambda_03_random.txt")

	t = np.linspace(0, len(data), len(data))

	plt.plot(t[0:-1:20], data[0:-1:20], ls='--', color='k', lw=2, label='Random selection')

	plt.xlabel(r"$t$")
	plt.ylabel(r"$I(t)$")

	plt.legend(loc='lower right')
	plt.show()

def prevalence_t_uniform_BA(ns):

	plt.figure(figsize=(8, 6))

	colors = ['g', 'r', 'b']
	markers = ['s', 'o', '^']

	i = -1

	for n in ns:

		i += 1

		data = np.loadtxt("Results/prevalence_BA_uniform_n_%s.txt" % str(n))

		t = np.linspace(0, len(data), len(data))

		#plt.plot(t, data, ls='-', color=colors[i], label=r'$n=%i$' % n)
		plt.plot(t[0:-1:20], data[0:-1:20], ls='-', color=colors[i], marker=markers[i], markersize=10, label=r'$n=%i$' % n)

	data = np.loadtxt("Results/prevalence_BA_random.txt")

	t = np.linspace(0, len(data), len(data))

	plt.plot(t[0:-1:20], data[0:-1:20], ls='--', color='k', lw=2, label='Random selection')

	plt.xlabel(r"$t$")
	plt.ylabel(r"$I(t)$")

	plt.legend(loc='lower right')
	plt.show()

def triplot_prevalence(ns):

	plt.subplots(figsize=(16, 8))

	colors = ['g', 'r', 'b']
	markers = ['s', 'o', '^']

	i = -1

	plt.subplot(1, 3, 1)
	for n in ns:

		i += 1

		data = np.loadtxt("Results/prevalence_BA_uniform_n_%s.txt" % str(n))

		t = np.linspace(0, len(data), len(data))

		#plt.plot(t, data, ls='-', color=colors[i], label=r'$n=%i$' % n)
		plt.plot(t[0:-1:20], data[0:-1:20], ls='-', color=colors[i], marker=markers[i], markersize=10, label=r'$n=%i$' % n)

	data = np.loadtxt("Results/prevalence_BA_random.txt")

	t = np.linspace(0, len(data), len(data))

	plt.plot(t[0:-1:20], data[0:-1:20], ls='--', color='k', lw=2, label='Random selection')

	plt.xlabel(r"$t$")
	plt.ylabel(r"$I(t)$")

	plt.legend(loc='lower right')
	
	plt.subplot(1, 3, 2)

	i = -1

	for n in ns:

		i += 1

		data = np.loadtxt("Results/small_world_lambda_03_n_%s.txt" % str(n))

		t = np.linspace(0, len(data), len(data))

		#plt.plot(t, data, ls='-', color=colors[i], label=r'$n=%i$' % n)
		plt.plot(t[0:-1:100], data[0:-1:100], ls='-', color=colors[i], marker=markers[i], markersize=10, label=r'$n=%i$' % n)

	data = np.loadtxt("Results/small_world_lambda_03_random.txt")

	t = np.linspace(0, len(data), len(data))

	plt.plot(t[0:-1:20], data[0:-1:20], ls='--', color='k', lw=2, label='Random selection')

	plt.xlabel(r"$t$")
	plt.ylabel(r"$I(t)$")

	plt.legend(loc='lower right')

	plt.subplot(1, 3, 3)

	i = -1

	for n in ns:

		i += 1

		data = np.loadtxt("Results/prevalence_WS_uniform_n_%s.txt" % str(n))

		t = np.linspace(0, len(data), len(data))

		#plt.plot(t, data, ls='-', color=colors[i], label=r'$n=%i$' % n)
		plt.plot(t[0:-1:20], data[0:-1:20], ls='-', color=colors[i], marker=markers[i], markersize=10, label=r'$n=%i$' % n)

	data = np.loadtxt("Results/prevalence_WS_random.txt")

	t = np.linspace(0, len(data), len(data))

	plt.plot(t[0:-1:20], data[0:-1:20], ls='--', color='k', lw=2, label='Random selection')

	plt.xlabel(r"$t$")
	plt.ylabel(r"$I(t)$")

	plt.legend(loc='lower right')
	plt.show()


def n_study(ks):

	markers = ['o', '^']

	i = -1

	plt.figure(figsize=(8, 6))

	for k in ks:

		i += 1

		data = np.loadtxt("Results/n_study_k_%s.txt" % k)

		prevalence = data[:, 0]
		errors = data[:, 1] / 2

		plt.ylim(0.32, 0.5)

		plt.errorbar(np.linspace(1, 10, 10), prevalence, ls='', marker=markers[i], markersize=9, markerfacecolor='None', yerr=errors, label=r'$k=%i$' % k)

	plt.xticks([i for i in range(1, 11)])

	plt.xlabel(r"$n$")
	plt.ylabel(r"$I$")

	plt.legend(loc="lower right")
	plt.show()

def gamma_study(ks):

	markers = ['o', '^']

	i = -1

	plt.figure(figsize=(8, 6))

	for k in ks:

		i += 1

		data = np.loadtxt("Results/gamma_study_k_%s.txt" % k)

		prevalence = data[:, 0]
		errors = data[:, 1] / 2

		plt.ylim(0.2, 0.5)

		plt.errorbar(np.arange(2, 4.1, 0.1), prevalence, ls='', marker=markers[i], markersize=9, markerfacecolor='None', yerr=errors, label=r'$k=%i$' % k)

	plt.xticks(np.arange(2, 4.2, 0.2))

	plt.xlabel(r"$n$")
	plt.ylabel(r"$I$")

	plt.legend(loc="lower right")
	plt.show()

def edge_density_study():
	
	data = np.loadtxt("Results/edge_density_study.txt")

	w = data[:, 0]
	E_SI = data[:, 1]
	E_SS = data[:, 2]
	E_II = data[:, 3]

	E_SI_upborder = 1 - E_II

	plt.plot(w, E_SS)
	plt.plot(w, E_SI_upborder)

	plt.fill_between(w, np.zeros(len(w)), E_SS)
	plt.fill_between(w, E_SS, E_SI_upborder)
	plt.fill_between(w, E_SI_upborder, np.ones(len(w)))

	xticks = np.arange(0, 2.2, 0.2)
	xlabels = []

	for i in range(len(xticks)):

		if i % 2 == 0:

			xlabels.append(str(round(xticks[i], 1)))

		else:

			xlabels.append('')

	yticks = np.arange(0, 1.1, 0.1)
	ylabels = []

	for i in range(len(yticks)):

		if i % 2 == 0:

			ylabels.append(str(round(yticks[i], 1)))

		else:

			ylabels.append('')

	plt.xticks(xticks, xlabels)
	plt.yticks(yticks, ylabels)

	plt.tick_params(labeltop=False, right=True)

	plt.xlim(0, 2)
	plt.ylim(0, 1)

	plt.show()


def difference_study():

	data = np.loadtxt("difference_study_RCM.txt")
	data_2 = np.loadtxt("difference_study_BA.txt")

	plt.plot(data[:, 0], data[:, 1], ls='', marker='o', label='RCM')
	plt.plot(data_2[:, 0]*2, data_2[:, 1], ls='', marker='^', label='BA') #BA <k> is double than input

	plt.legend()
	plt.show()


names = ["Single", "Uniform", "Poisson"]

#alpha_study_SIS_WS(names)
#alpha_study_SIR_WS(names)
#alpha_study_SIR_BA(names)

#lambda_study_SIS_WS(names)
#lambda_study_SIR_WS(names)

lambdas = [0.024, 0.025, 0.028]

#dR_lambda(lambdas)


################
#	 PAPER 2   #
################

ns = [1, 5, 10]

gammas = [2, 3, 4]

ks = [5, 10]

#prevalence_t_uniform(ns)
#prevalence_t_powerlaw(gammas)

#prevalence_t_uniform_BA(ns)
#prevalence_t_uniform_WS(ns)
#prevalence_t_small_world(ns)

#triplot_prevalence(ns)


#n_study(ks)
#gamma_study(ks)
#edge_density_study()

difference_study()





