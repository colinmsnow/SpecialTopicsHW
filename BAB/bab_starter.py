import picos as pic
from picos import RealVariable
from copy import deepcopy
from heapq import *
import heapq as hq
import numpy as np
import itertools
import math
counter = itertools.count() 

class BBTreeNode():
    def __init__(self, vars = [], constraints = [], objective='', prob=None):
        self.vars = vars
        self.constraints = constraints
        self.objective = objective
        self.prob = prob

    def __deepcopy__(self, memo):
        '''
        Deepcopies the picos problem
        This overrides the system's deepcopy method bc it doesn't work on classes by itself
        '''
        newprob = pic.Problem.clone(self.prob)
        return BBTreeNode(self.vars, newprob.constraints, self.objective, newprob)
    
    def buildProblem(self):
        '''
        Bulids the initial Picos problem
        '''
        prob=pic.Problem()

        # print(self.constraints)
        prob.add_list_of_constraints(self.constraints)    
        
        prob.set_objective('max', self.objective)
        self.prob = prob
        return self.prob

    def is_integral(self):
        '''
        Checks if all variables (excluding the one we're maxing) are integers
        '''
        for v in self.vars[:-1]:
            if v.value == None or abs(round(v.value) - float(v.value)) > 1e-4 :
                return False
        return True

    def branch_floor(self, branch_var):
        '''
        Makes a child where xi <= floor(xi)
        '''
        n1 = deepcopy(self)
        n1.prob.add_constraint( branch_var <= math.floor(branch_var.value) ) # add in the new binary constraint

        return n1

    def branch_ceil(self, branch_var):
        '''
        Makes a child where xi >= ceiling(xi)
        '''
        n2 = deepcopy(self)
        n2.prob.add_constraint( branch_var >= math.ceil(branch_var.value) ) # add in the new binary constraint
        return n2



    def first_noninteger(self):
        for v in self.vars[:-1]:
            if v.value == None or abs(round(v.value) - float(v.value)) > 1e-4 :
                return v
        return None


    def bbsolve(self):
        '''
        Use the branch and bound method to solve an integer program
        This function should return:
            return bestres, bestnode_vars

        where bestres = value of the maximized objective function
              bestnode_vars = the list of variables that create bestres
        '''


        # print(self)
        # print(self.vars)
        # print(self.constraints)
        # print(self.prob)
        # quit()

        # these lines build up the initial problem and adds it to a heap
        root = self
        res = root.buildProblem().solve(solver='cvxopt')
        heap = [(res, next(counter), root)]
        bestres = -1e20 # a small arbitrary initial best objective value
        bestnode_vars = root.vars # initialize bestnode_vars to the root vars

        if self.is_integral():
            return res.value, bestnode_vars

        bestres, bestnode_vars = recursive_solve(self, bestres, bestnode_vars)




        # quit()
        # print(bestres, bestnode_vars)
        return bestres, bestnode_vars

def recursive_solve(obj, bestres, bestnode_vars):

    # soln = obj.buildProblem().solve(solver='cvxopt')
    try:
        soln = obj.prob.solve(solver='cvxopt')
        # print(soln.value)

    except pic.modeling.problem.SolutionFailure:
        # no valid solution
        # print("FAILED SOLUTION")
        return 0, bestnode_vars

    if obj.is_integral(): # If everything is integer we are done
        # print("IS INTEGRAL")

        bestres = round(soln.value)
        bestnode_vars = [round(i) for i in obj.vars]
        # print(bestres, bestnode_vars)
        return bestres, bestnode_vars

    if soln.value > bestres: # don't bother if the solution is worse

        # print(soln.value)

        branch_var = obj.first_noninteger()

        res, node_vars = recursive_solve(obj.branch_floor(branch_var), bestres, bestnode_vars)
        res2, node_vars2 = recursive_solve(obj.branch_ceil(branch_var), bestres, bestnode_vars)

        if res > bestres:
            bestres = res
            bestnode_vars = node_vars

        if res2 > bestres:
            bestres = res2
            bestnode_vars = node_vars2

    return bestres, bestnode_vars


 
