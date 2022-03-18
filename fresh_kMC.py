import numpy
import matplotlib.pyplot as plt
import gillespy2
from gillespy2 import Model, Species, Parameter, Reaction
from gillespy2.solvers.numpy.basic_ode_solver import BasicODESolver

class ToggleSwitch(Model):
    def __init__(self, parameter_values = None):
        # Initialize the model.
        Model.__init__(self, name = "toggle_switch")

        # Define parameters.
        alpha1 = Parameter(name = 'alpha1', expression = 1)
        alpha2 = Parameter(name = 'alpha2', expression = 1)
        beta   = Parameter(name = 'beta',   expression = 2.0)
        gamma  = Parameter(name = 'gamma',  expression = 2.0)
        mu     = Parameter(name = 'mu',     expression = 1.0)
        self.add_parameter([alpha1, alpha2, beta, gamma, mu])

        # Define molecular species.
        U = Species(name = 'U', initial_value = 10)
        V = Species(name = 'V', initial_value = 10)
        self.add_species([U, V])

        # Define reactions.
        cu = Reaction(name = "r1", reactants = {}, products = {U:1},
                      propensity_function = "alpha1/(1+pow(V,beta))")
        cv = Reaction(name = "r2", reactants = {}, products = {V:1},
                      propensity_function = "alpha2/(1+pow(U,gamma))")
        du = Reaction(name = "r3", reactants = {U:1}, products = {},
                      rate = mu)
        dv = Reaction(name = "r4", reactants = {V:1}, products = {},
                      rate = mu)
        self.add_reaction([cu, cv, du, dv])
        self.timespan(numpy.linspace(0, 100, 101))
        
model = ToggleSwitch()
s_results = model.run(show_labels = False)
d_results = model.run(solver = BasicODESolver, show_labels = False)
plt.plot(s_results[0][:,0], s_results[0][:,1], '-r', label='U')
plt.plot(s_results[0][:,0], s_results[0][:,2], '-b', label='V')
plt.plot([0], [11])
plt.title('Stochastic Switch')
plt.legend(loc = 'best')
plt.show()