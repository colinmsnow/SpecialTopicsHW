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
        """ Gets the first variable in the eq that is not an integer
            It could be any noninteger but the first is easy"""
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

        # these lines build up the initial problem and adds it to a heap
        root = self
        res = root.buildProblem().solve(solver='cvxopt')
        bestres = -1e20 # a small arbitrary initial best objective value
        bestnode_vars = root.vars # initialize bestnode_vars to the root vars

        bestres, bestnode_vars = recursive_solve(self, bestres, bestnode_vars)

        return bestres, bestnode_vars


def recursive_solve(obj, bestres, bestnode_vars):
    """ Recursively iterate through possible solutions, adding
        constraints along the way until all variables are integer.
        Return max all-integer solution """

    try:
        soln = obj.prob.solve(solver='cvxopt')

    except pic.modeling.problem.SolutionFailure:
        # no valid solution
        return 0, bestnode_vars

    if obj.is_integral(): # If everything is integer we are done

        res = round(soln.value)
        node_vars = [round(i) for i in obj.vars]
        return res, node_vars

    if soln.value > bestres: # don't bother if the solution is worse

        branch_var = obj.first_noninteger() # Get a variable to branch on

        # Take the upper and lower branch
        res, node_vars = recursive_solve(obj.branch_floor(branch_var), bestres, bestnode_vars)
        res2, node_vars2 = recursive_solve(obj.branch_ceil(branch_var), bestres, bestnode_vars)

        # update best solution if necessary
        if res > bestres:
            bestres = res
            bestnode_vars = node_vars

        if res2 > bestres:
            bestres = res2
            bestnode_vars = node_vars2

    return bestres, bestnode_vars


 
