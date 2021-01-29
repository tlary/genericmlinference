from opossum import UserInterface
import pandas as pd

# number of observations N and number of covariates k
N = 20000
k = 20
# initilizing class
# 10 binary variables
u = UserInterface(N, k, seed=None, categorical_covariates = [10,2])

# assign treatment and generate treatment effect inside of class object
u.generate_treatment(random_assignment = True,
                     assignment_prob = 0.5,
                     constant_pos = False,
                     constant_neg = False,
                     heterogeneous_pos = True,
                     heterogeneous_neg = True,
                     no_treatment = False,
                     discrete_heterogeneous = False,
                     treatment_option_weights = [0, 0, 0.5, 0.5, 0, 0],
                     intensity = 10)
# generate output variable y and return all 4 variables
y, X, assignment, treatment = u.output_data(binary=False, x_y_relation = 'linear_simple')

pd.DataFrame(y).to_csv("y.csv")
pd.DataFrame(X).to_csv("X.csv")
pd.DataFrame(assignment).to_csv("assignment.csv")
pd.DataFrame(treatment).to_csv("treatment.csv")

