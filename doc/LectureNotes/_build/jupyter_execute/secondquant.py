#!/usr/bin/env python
# coding: utf-8

# <!-- HTML file automatically generated from DocOnce source (https://github.com/doconce/doconce/)
# doconce format html secondquant.do.txt  -->

# # Introduction to  second quantization
# 
# We introduce the time-independent  operators
# $a_\alpha^{\dagger}$ and $a_\alpha$   which create and annihilate, respectively, a particle 
# in the single-particle state 
# $\varphi_\alpha$. 
# We define the fermion creation operator
# $a_\alpha^{\dagger}$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-1a"></div>
# 
# $$
# \begin{equation}
# 	a_\alpha^{\dagger}|0\rangle \equiv  |\alpha\rangle  \label{eq:2-1a} \tag{1},
# \end{equation}
# $$

# and

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-1b"></div>
# 
# $$
# \begin{equation}
# 	a_\alpha^{\dagger}|\alpha_1\dots \alpha_n\rangle_{\mathrm{AS}} \equiv  |\alpha\alpha_1\dots \alpha_n\rangle_{\mathrm{AS}} \label{eq:2-1b} \tag{2}
# \end{equation}
# $$

# In Eq. ([1](#eq:2-1a)) 
# the operator  $a_\alpha^{\dagger}$  acts on the vacuum state 
# $|0\rangle$, which does not contain any particles. Alternatively, we could define  a closed-shell nucleus or atom as our new vacuum, but then
# we need to introduce the particle-hole  formalism, see the discussion to come. 
# 
# In Eq. ([2](#eq:2-1b)) $a_\alpha^{\dagger}$ acts on an antisymmetric $n$-particle state and 
# creates an antisymmetric $(n+1)$-particle state, where the one-body state 
# $\varphi_\alpha$ is occupied, under the condition that
# $\alpha \ne \alpha_1, \alpha_2, \dots, \alpha_n$. 
# It follows that we can express an antisymmetric state as the product of the creation
# operators acting on the vacuum state.

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-2"></div>
# 
# $$
# \begin{equation}
# 	|\alpha_1\dots \alpha_n\rangle_{\mathrm{AS}} = a_{\alpha_1}^{\dagger} a_{\alpha_2}^{\dagger} \dots a_{\alpha_n}^{\dagger} |0\rangle \label{eq:2-2} \tag{3}
# \end{equation}
# $$

# It is easy to derive the commutation and anticommutation rules  for the fermionic creation operators 
# $a_\alpha^{\dagger}$. Using the antisymmetry of the states 
# ([3](#eq:2-2))

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-3a"></div>
# 
# $$
# \begin{equation}
# 	|\alpha_1\dots \alpha_i\dots \alpha_k\dots \alpha_n\rangle_{\mathrm{AS}} = 
# 		- |\alpha_1\dots \alpha_k\dots \alpha_i\dots \alpha_n\rangle_{\mathrm{AS}} \label{eq:2-3a} \tag{4}
# \end{equation}
# $$

# we obtain

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-3b"></div>
# 
# $$
# \begin{equation}
# 	 a_{\alpha_i}^{\dagger}  a_{\alpha_k}^{\dagger} = - a_{\alpha_k}^{\dagger} a_{\alpha_i}^{\dagger} \label{eq:2-3b} \tag{5}
# \end{equation}
# $$

# Using the Pauli principle

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-4a"></div>
# 
# $$
# \begin{equation}
# 	|\alpha_1\dots \alpha_i\dots \alpha_i\dots \alpha_n\rangle_{\mathrm{AS}} = 0 \label{eq:2-4a} \tag{6}
# \end{equation}
# $$

# it follows that

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-4b"></div>
# 
# $$
# \begin{equation}
# 	a_{\alpha_i}^{\dagger}  a_{\alpha_i}^{\dagger} = 0. \label{eq:2-4b} \tag{7}
# \end{equation}
# $$

# If we combine Eqs. ([5](#eq:2-3b)) and ([7](#eq:2-4b)), we obtain the well-known anti-commutation rule

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-5"></div>
# 
# $$
# \begin{equation}
# 	a_{\alpha}^{\dagger}  a_{\beta}^{\dagger} + a_{\beta}^{\dagger}  a_{\alpha}^{\dagger} \equiv 
# 		\{a_{\alpha}^{\dagger},a_{\beta}^{\dagger}\} = 0 \label{eq:2-5} \tag{8}
# \end{equation}
# $$

# The hermitian conjugate  of $a_\alpha^{\dagger}$ is

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-6"></div>
# 
# $$
# \begin{equation}
# 	a_{\alpha} = ( a_{\alpha}^{\dagger} )^{\dagger} \label{eq:2-6} \tag{9}
# \end{equation}
# $$

# If we take the hermitian conjugate of Eq. ([8](#eq:2-5)), we arrive at

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-7"></div>
# 
# $$
# \begin{equation}
# 	\{a_{\alpha},a_{\beta}\} = 0 \label{eq:2-7} \tag{10}
# \end{equation}
# $$

# What is the physical interpretation of the operator $a_\alpha$ and what is the effect of 
# $a_\alpha$ on a given state $|\alpha_1\alpha_2\dots\alpha_n\rangle_{\mathrm{AS}}$? 
# Consider the following matrix element

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-8"></div>
# 
# $$
# \begin{equation}
# 	\langle\alpha_1\alpha_2 \dots \alpha_n|a_\alpha|\alpha_1'\alpha_2' \dots \alpha_m'\rangle \label{eq:2-8} \tag{11}
# \end{equation}
# $$

# where both sides are antisymmetric. We  distinguish between two cases. The first (1) is when
# $\alpha \in \{\alpha_i\}$. Using the Pauli principle of Eq. ([6](#eq:2-4a)) it follows

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-9a"></div>
# 
# $$
# \begin{equation}
# 		\langle\alpha_1\alpha_2 \dots \alpha_n|a_\alpha = 0 \label{eq:2-9a} \tag{12}
# \end{equation}
# $$

# The second (2) case is when $\alpha \notin \{\alpha_i\}$. It follows that an hermitian conjugation

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-9b"></div>
# 
# $$
# \begin{equation}
# 		\langle \alpha_1\alpha_2 \dots \alpha_n|a_\alpha = \langle\alpha\alpha_1\alpha_2 \dots \alpha_n|  \label{eq:2-9b} \tag{13}
# \end{equation}
# $$

# Eq. ([13](#eq:2-9b)) holds for case (1) since the lefthand side is zero due to the Pauli principle. We write
# Eq. ([11](#eq:2-8)) as

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-10"></div>
# 
# $$
# \begin{equation}
# 	\langle\alpha_1\alpha_2 \dots \alpha_n|a_\alpha|\alpha_1'\alpha_2' \dots \alpha_m'\rangle = 
# 	\langle \alpha_1\alpha_2 \dots \alpha_n|\alpha\alpha_1'\alpha_2' \dots \alpha_m'\rangle \label{eq:2-10} \tag{14}
# \end{equation}
# $$

# Here we must have $m = n+1$ if Eq. ([14](#eq:2-10)) has to be trivially different from zero.
# 
# For the last case, the minus and plus signs apply when the sequence 
# $\alpha ,\alpha_1, \alpha_2, \dots, \alpha_n$ and 
# $\alpha_1', \alpha_2', \dots, \alpha_{n+1}'$ are related to each other via even and odd permutations.
# If we assume that  $\alpha \notin \{\alpha_i\}$ we obtain

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-12"></div>
# 
# $$
# \begin{equation}
# 	\langle\alpha_1\alpha_2 \dots \alpha_n|a_\alpha|\alpha_1'\alpha_2' \dots \alpha_{n+1}'\rangle = 0 \label{eq:2-12} \tag{15}
# \end{equation}
# $$

# when $\alpha \in \{\alpha_i'\}$. If $\alpha \notin \{\alpha_i'\}$, we obtain

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-13a"></div>
# 
# $$
# \begin{equation}
# 	a_\alpha\underbrace{|\alpha_1'\alpha_2' \dots \alpha_{n+1}'}\rangle_{\neq \alpha} = 0 \label{eq:2-13a} \tag{16}
# \end{equation}
# $$

# and in particular

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-13b"></div>
# 
# $$
# \begin{equation}
# 	a_\alpha |0\rangle = 0 \label{eq:2-13b} \tag{17}
# \end{equation}
# $$

# If $\{\alpha\alpha_i\} = \{\alpha_i'\}$, performing the right permutations, the sequence
# $\alpha ,\alpha_1, \alpha_2, \dots, \alpha_n$ is identical with the sequence
# $\alpha_1', \alpha_2', \dots, \alpha_{n+1}'$. This results in

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-14"></div>
# 
# $$
# \begin{equation}
# 	\langle\alpha_1\alpha_2 \dots \alpha_n|a_\alpha|\alpha\alpha_1\alpha_2 \dots \alpha_{n}\rangle = 1 \label{eq:2-14} \tag{18}
# \end{equation}
# $$

# and thus

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-15"></div>
# 
# $$
# \begin{equation}
# 	a_\alpha |\alpha\alpha_1\alpha_2 \dots \alpha_{n}\rangle = |\alpha_1\alpha_2 \dots \alpha_{n}\rangle \label{eq:2-15} \tag{19}
# \end{equation}
# $$

# The action of the operator 
# $a_\alpha$ from the left on a state vector  is to to remove  one particle in the state
# $\alpha$. 
# If the state vector does not contain the single-particle state $\alpha$, the outcome of the operation is zero.
# The operator  $a_\alpha$ is normally called for a destruction or annihilation operator.
# 
# The next step is to establish the  commutator algebra of $a_\alpha^{\dagger}$ and
# $a_\beta$. 
# 
# The action of the anti-commutator 
# $\{a_\alpha^{\dagger}$,$a_\alpha\}$ on a given $n$-particle state is

# $$
# a_\alpha^{\dagger} a_\alpha \underbrace{|\alpha_1\alpha_2 \dots \alpha_{n}\rangle}_{\neq \alpha} = 0 \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-16a"></div>
# 
# $$
# \begin{equation} 
# 	a_\alpha a_\alpha^{\dagger} \underbrace{|\alpha_1\alpha_2 \dots \alpha_{n}\rangle}_{\neq \alpha} =
# 	a_\alpha \underbrace{|\alpha \alpha_1\alpha_2 \dots \alpha_{n}\rangle}_{\neq \alpha} = 
# 	\underbrace{|\alpha_1\alpha_2 \dots \alpha_{n}\rangle}_{\neq \alpha} \label{eq:2-16a} \tag{20}
# \end{equation}
# $$

# if the single-particle state $\alpha$ is not contained in the state.
# 
#  If it is present
# we arrive at

# $$
# a_\alpha^{\dagger} a_\alpha |\alpha_1\alpha_2 \dots \alpha_{k}\alpha \alpha_{k+1} \dots \alpha_{n-1}\rangle =
# 	a_\alpha^{\dagger} a_\alpha (-1)^k |\alpha \alpha_1\alpha_2 \dots \alpha_{n-1}\rangle \nonumber
# $$

# $$
# = (-1)^k |\alpha \alpha_1\alpha_2 \dots \alpha_{n-1}\rangle =
# 	|\alpha_1\alpha_2 \dots \alpha_{k}\alpha \alpha_{k+1} \dots \alpha_{n-1}\rangle \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-16b"></div>
# 
# $$
# \begin{equation} 
# 	a_\alpha a_\alpha^{\dagger}|\alpha_1\alpha_2 \dots \alpha_{k}\alpha \alpha_{k+1} \dots \alpha_{n-1}\rangle = 0 \label{eq:2-16b} \tag{21}
# \end{equation}
# $$

# From Eqs. ([20](#eq:2-16a)) and  ([21](#eq:2-16b)) we arrive at

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-17"></div>
# 
# $$
# \begin{equation}
# 	\{a_\alpha^{\dagger} , a_\alpha \} = a_\alpha^{\dagger} a_\alpha + a_\alpha a_\alpha^{\dagger} = 1 \label{eq:2-17} \tag{22}
# \end{equation}
# $$

# The action of $\left\{a_\alpha^{\dagger}, a_\beta\right\}$, with 
# $\alpha \ne \beta$ on a given state yields three possibilities. 
# The first case is a state vector which contains both $\alpha$ and $\beta$, then either 
# $\alpha$ or $\beta$ and finally none of them.
# 
# The first case results in

# $$
# a_\alpha^{\dagger} a_\beta |\alpha\beta\alpha_1\alpha_2 \dots \alpha_{n-2}\rangle = 0 \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-18a"></div>
# 
# $$
# \begin{equation} 
# 	a_\beta a_\alpha^{\dagger} |\alpha\beta\alpha_1\alpha_2 \dots \alpha_{n-2}\rangle = 0 \label{eq:2-18a} \tag{23}
# \end{equation}
# $$

# while the second case gives

# $$
# a_\alpha^{\dagger} a_\beta |\beta \underbrace{\alpha_1\alpha_2 \dots \alpha_{n-1}}_{\neq \alpha}\rangle = 
# 	 	|\alpha \underbrace{\alpha_1\alpha_2 \dots \alpha_{n-1}}_{\neq  \alpha}\rangle \nonumber
# $$

# $$
# a_\beta a_\alpha^{\dagger} |\beta \underbrace{\alpha_1\alpha_2 \dots \alpha_{n-1}}_{\neq \alpha}\rangle =
# 		a_\beta |\alpha\beta\underbrace{\beta \alpha_1\alpha_2 \dots \alpha_{n-1}}_{\neq \alpha}\rangle \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-18b"></div>
# 
# $$
# \begin{equation} 
# 	= - |\alpha\underbrace{\alpha_1\alpha_2 \dots \alpha_{n-1}}_{\neq \alpha}\rangle \label{eq:2-18b} \tag{24}
# \end{equation}
# $$

# Finally if the state vector does not contain $\alpha$ and $\beta$

# $$
# a_\alpha^{\dagger} a_\beta |\underbrace{\alpha_1\alpha_2 \dots \alpha_{n}}_{\neq \alpha,\beta}\rangle = 0 \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-18c"></div>
# 
# $$
# \begin{equation} 
# 	a_\beta a_\alpha^{\dagger} |\underbrace{\alpha_1\alpha_2 \dots \alpha_{n}}_{\neq \alpha,\beta}\rangle = 
# 		a_\beta |\alpha \underbrace{\alpha_1\alpha_2 \dots \alpha_{n}}_{\neq \alpha,\beta}\rangle = 0 \label{eq:2-18c} \tag{25}
# \end{equation}
# $$

# For all three cases we have

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-19"></div>
# 
# $$
# \begin{equation}
# 	\{a_\alpha^{\dagger},a_\beta \} = a_\alpha^{\dagger} a_\beta + a_\beta a_\alpha^{\dagger} = 0, \quad \alpha \neq \beta \label{eq:2-19} \tag{26}
# \end{equation}
# $$

# We can summarize  our findings in Eqs. ([22](#eq:2-17)) and ([26](#eq:2-19)) as

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-20"></div>
# 
# $$
# \begin{equation}
# 	\{a_\alpha^{\dagger},a_\beta \} = \delta_{\alpha\beta} \label{eq:2-20} \tag{27}
# \end{equation}
# $$

# with  $\delta_{\alpha\beta}$ is the Kroenecker $\delta$-symbol.
# 
# The properties of the creation and annihilation operators can be summarized as (for fermions)

# $$
# a_\alpha^{\dagger}|0\rangle \equiv  |\alpha\rangle,
# $$

# and

# $$
# a_\alpha^{\dagger}|\alpha_1\dots \alpha_n\rangle_{\mathrm{AS}} \equiv  |\alpha\alpha_1\dots \alpha_n\rangle_{\mathrm{AS}}.
# $$

# from which follows

# $$
# |\alpha_1\dots \alpha_n\rangle_{\mathrm{AS}} = a_{\alpha_1}^{\dagger} a_{\alpha_2}^{\dagger} \dots a_{\alpha_n}^{\dagger} |0\rangle.
# $$

# The hermitian conjugate has the folowing properties

# $$
# a_{\alpha} = ( a_{\alpha}^{\dagger} )^{\dagger}.
# $$

# Finally we found

# $$
# a_\alpha\underbrace{|\alpha_1'\alpha_2' \dots \alpha_{n+1}'}\rangle_{\neq \alpha} = 0, \quad
# 		\textrm{in particular } a_\alpha |0\rangle = 0,
# $$

# and

# $$
# a_\alpha |\alpha\alpha_1\alpha_2 \dots \alpha_{n}\rangle = |\alpha_1\alpha_2 \dots \alpha_{n}\rangle,
# $$

# and the corresponding commutator algebra

# $$
# \{a_{\alpha}^{\dagger},a_{\beta}^{\dagger}\} = \{a_{\alpha},a_{\beta}\} = 0 \hspace{0.5cm} 
# \{a_\alpha^{\dagger},a_\beta \} = \delta_{\alpha\beta}.
# $$

# ## One-body operators in second quantization
# 
# A very useful operator is the so-called number-operator.  Most physics cases  we will
# study in this text conserve the total number of particles.  The number operator is therefore
# a useful quantity which allows us to test that our many-body formalism  conserves the number of particles.
# In for example $(d,p)$ or $(p,d)$ reactions it is important to be able to describe quantum mechanical states
# where particles get added or removed.
# A creation operator $a_\alpha^{\dagger}$ adds one particle to the single-particle state
# $\alpha$ of a give many-body state vector, while an annihilation operator $a_\alpha$ 
# removes a particle from a single-particle
# state $\alpha$. 
# 
# Let us consider an operator proportional with $a_\alpha^{\dagger} a_\beta$ and 
# $\alpha=\beta$. It acts on an $n$-particle state 
# resulting in

# <!-- Equation labels as ordinary links -->
# <div id="_auto1"></div>
# 
# $$
# \begin{equation}
# 	a_\alpha^{\dagger} a_\alpha |\alpha_1\alpha_2 \dots \alpha_{n}\rangle = 
# 	\begin{cases}
# 		0  &\alpha \notin \{\alpha_i\} \\
# 		\\
# 		|\alpha_1\alpha_2 \dots \alpha_{n}\rangle & \alpha \in \{\alpha_i\}
# 	\end{cases}
# \label{_auto1} \tag{28}
# \end{equation}
# $$

# Summing over all possible one-particle states we arrive at

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-21"></div>
# 
# $$
# \begin{equation}
# 	\left( \sum_\alpha a_\alpha^{\dagger} a_\alpha \right) |\alpha_1\alpha_2 \dots \alpha_{n}\rangle = 
# 	n |\alpha_1\alpha_2 \dots \alpha_{n}\rangle \label{eq:2-21} \tag{29}
# \end{equation}
# $$

# The operator

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-22"></div>
# 
# $$
# \begin{equation}
# 	\hat{N} = \sum_\alpha a_\alpha^{\dagger} a_\alpha \label{eq:2-22} \tag{30}
# \end{equation}
# $$

# is called the number operator since it counts the number of particles in a give state vector when it acts 
# on the different single-particle states.  It acts on one single-particle state at the time and falls 
# therefore under category one-body operators.
# Next we look at another important one-body operator, namely $\hat{H}_0$ and study its operator form in the 
# occupation number representation.
# 
# We want to obtain an expression for a one-body operator which conserves the number of particles.
# Here we study the one-body operator for the kinetic energy plus an eventual external one-body potential.
# The action of this operator on a particular $n$-body state with its pertinent expectation value has already
# been studied in coordinate  space.
# In coordinate space the operator reads

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-23"></div>
# 
# $$
# \begin{equation}
# 	\hat{H}_0 = \sum_i \hat{h}_0(x_i) \label{eq:2-23} \tag{31}
# \end{equation}
# $$

# and the anti-symmetric $n$-particle Slater determinant is defined as

# $$
# \Phi(x_1, x_2,\dots ,x_n,\alpha_1,\alpha_2,\dots, \alpha_n)= \frac{1}{\sqrt{n!}} \sum_p (-1)^p\hat{P}\psi_{\alpha_1}(x_1)\psi_{\alpha_2}(x_2) \dots \psi_{\alpha_n}(x_n).
# $$

# Defining

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-25"></div>
# 
# $$
# \begin{equation}
# 	\hat{h}_0(x_i) \psi_{\alpha_i}(x_i) = \sum_{\alpha_k'} \psi_{\alpha_k'}(x_i) \langle\alpha_k'|\hat{h}_0|\alpha_k\rangle \label{eq:2-25} \tag{32}
# \end{equation}
# $$

# we can easily  evaluate the action of $\hat{H}_0$ on each product of one-particle functions in Slater determinant.
# From Eq. ([32](#eq:2-25))  we obtain the following result without  permuting any particle pair

# $$
# \left( \sum_i \hat{h}_0(x_i) \right) \psi_{\alpha_1}(x_1)\psi_{\alpha_2}(x_2) \dots \psi_{\alpha_n}(x_n) \nonumber
# $$

# $$
# =\sum_{\alpha_1'} \langle \alpha_1'|\hat{h}_0|\alpha_1\rangle 
# 		\psi_{\alpha_1'}(x_1)\psi_{\alpha_2}(x_2) \dots \psi_{\alpha_n}(x_n) \nonumber
# $$

# $$
# +\sum_{\alpha_2'} \langle \alpha_2'|\hat{h}_0|\alpha_2\rangle
# 		\psi_{\alpha_1}(x_1)\psi_{\alpha_2'}(x_2) \dots \psi_{\alpha_n}(x_n) \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-26"></div>
# 
# $$
# \begin{equation} 
# 	+\sum_{\alpha_n'} \langle \alpha_n'|\hat{h}_0|\alpha_n\rangle
# 		\psi_{\alpha_1}(x_1)\psi_{\alpha_2}(x_2) \dots \psi_{\alpha_n'}(x_n) \label{eq:2-26} \tag{33}
# \end{equation}
# $$

# If we interchange particles $1$ and $2$  we obtain

# $$
# \left( \sum_i \hat{h}_0(x_i) \right) \psi_{\alpha_1}(x_2)\psi_{\alpha_1}(x_2) \dots \psi_{\alpha_n}(x_n) \nonumber
# $$

# $$
# =\sum_{\alpha_2'} \langle \alpha_2'|\hat{h}_0|\alpha_2\rangle 
# 		\psi_{\alpha_1}(x_2)\psi_{\alpha_2'}(x_1) \dots \psi_{\alpha_n}(x_n) \nonumber
# $$

# $$
# +\sum_{\alpha_1'} \langle \alpha_1'|\hat{h}_0|\alpha_1\rangle
# 		\psi_{\alpha_1'}(x_2)\psi_{\alpha_2}(x_1) \dots \psi_{\alpha_n}(x_n) \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-27"></div>
# 
# $$
# \begin{equation} 
# 	+\sum_{\alpha_n'} \langle \alpha_n'|\hat{h}_0|\alpha_n\rangle
# 		\psi_{\alpha_1}(x_2)\psi_{\alpha_1}(x_2) \dots \psi_{\alpha_n'}(x_n) \label{eq:2-27} \tag{34}
# \end{equation}
# $$

# We can continue by computing all possible permutations. 
# We rewrite also our Slater determinant in its second quantized form and skip the dependence on the quantum numbers $x_i.$
# Summing up all contributions and taking care of all phases
# $(-1)^p$ we arrive at

# $$
# \hat{H}_0|\alpha_1,\alpha_2,\dots, \alpha_n\rangle = \sum_{\alpha_1'}\langle \alpha_1'|\hat{h}_0|\alpha_1\rangle
# 		|\alpha_1'\alpha_2 \dots \alpha_{n}\rangle \nonumber
# $$

# $$
# + \sum_{\alpha_2'} \langle \alpha_2'|\hat{h}_0|\alpha_2\rangle
# 		|\alpha_1\alpha_2' \dots \alpha_{n}\rangle \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-28"></div>
# 
# $$
# \begin{equation} 
# 	+ \sum_{\alpha_n'} \langle \alpha_n'|\hat{h}_0|\alpha_n\rangle
# 		|\alpha_1\alpha_2 \dots \alpha_{n}'\rangle \label{eq:2-28} \tag{35}
# \end{equation}
# $$

# In Eq. ([35](#eq:2-28)) 
# we have expressed the action of the one-body operator
# of Eq. ([31](#eq:2-23)) on the  $n$-body state in its second quantized form.
# This equation can be further manipulated if we use the properties of the creation and annihilation operator
# on each primed quantum number, that is

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-29"></div>
# 
# $$
# \begin{equation}
# 	|\alpha_1\alpha_2 \dots \alpha_k' \dots \alpha_{n}\rangle = 
# 		a_{\alpha_k'}^{\dagger}  a_{\alpha_k} |\alpha_1\alpha_2 \dots \alpha_k \dots \alpha_{n}\rangle \label{eq:2-29} \tag{36}
# \end{equation}
# $$

# Inserting this in the right-hand side of Eq. ([35](#eq:2-28)) results in

# $$
# \hat{H}_0|\alpha_1\alpha_2 \dots \alpha_{n}\rangle = \sum_{\alpha_1'}\langle \alpha_1'|\hat{h}_0|\alpha_1\rangle
# 		a_{\alpha_1'}^{\dagger}  a_{\alpha_1} |\alpha_1\alpha_2 \dots \alpha_{n}\rangle \nonumber
# $$

# $$
# + \sum_{\alpha_2'} \langle \alpha_2'|\hat{h}_0|\alpha_2\rangle
# 		a_{\alpha_2'}^{\dagger}  a_{\alpha_2} |\alpha_1\alpha_2 \dots \alpha_{n}\rangle \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# $$
# + \sum_{\alpha_n'} \langle \alpha_n'|\hat{h}_0|\alpha_n\rangle
# 		a_{\alpha_n'}^{\dagger}  a_{\alpha_n} |\alpha_1\alpha_2 \dots \alpha_{n}\rangle \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-30a"></div>
# 
# $$
# \begin{equation} 
# 	= \sum_{\alpha, \beta} \langle \alpha|\hat{h}_0|\beta\rangle a_\alpha^{\dagger} a_\beta 
# 		|\alpha_1\alpha_2 \dots \alpha_{n}\rangle \label{eq:2-30a} \tag{37}
# \end{equation}
# $$

# In the number occupation representation or second quantization we get the following expression for a one-body 
# operator which conserves the number of particles

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-30b"></div>
# 
# $$
# \begin{equation}
# 	\hat{H}_0 = \sum_{\alpha\beta} \langle \alpha|\hat{h}_0|\beta\rangle a_\alpha^{\dagger} a_\beta \label{eq:2-30b} \tag{38}
# \end{equation}
# $$

# Obviously, $\hat{H}_0$ can be replaced by any other one-body  operator which preserved the number
# of particles. The stucture of the operator is therefore not limited to say the kinetic or single-particle energy only.
# 
# The opearator $\hat{H}_0$ takes a particle from the single-particle state $\beta$  to the single-particle state $\alpha$ 
# with a probability for the transition given by the expectation value $\langle \alpha|\hat{h}_0|\beta\rangle$.
# 
# It is instructive to verify Eq. ([38](#eq:2-30b)) by computing the expectation value of $\hat{H}_0$ 
# between two single-particle states

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-30c"></div>
# 
# $$
# \begin{equation}
# 	\langle \alpha_1|\hat{h}_0|\alpha_2\rangle = \sum_{\alpha\beta} \langle \alpha|\hat{h}_0|\beta\rangle
# 		\langle 0|a_{\alpha_1}a_\alpha^{\dagger} a_\beta a_{\alpha_2}^{\dagger}|0\rangle \label{eq:2-30c} \tag{39}
# \end{equation}
# $$

# Using the commutation relations for the creation and annihilation operators we have

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-30d"></div>
# 
# $$
# \begin{equation}
# a_{\alpha_1}a_\alpha^{\dagger} a_\beta a_{\alpha_2}^{\dagger} = (\delta_{\alpha \alpha_1} - a_\alpha^{\dagger} a_{\alpha_1} )(\delta_{\beta \alpha_2} - a_{\alpha_2}^{\dagger} a_{\beta} ), \label{eq:2-30d} \tag{40}
# \end{equation}
# $$

# which results in

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-30e"></div>
# 
# $$
# \begin{equation}
# \langle 0|a_{\alpha_1}a_\alpha^{\dagger} a_\beta a_{\alpha_2}^{\dagger}|0\rangle = \delta_{\alpha \alpha_1} \delta_{\beta \alpha_2} \label{eq:2-30e} \tag{41}
# \end{equation}
# $$

# and

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-30f"></div>
# 
# $$
# \begin{equation}
# \langle \alpha_1|\hat{h}_0|\alpha_2\rangle = \sum_{\alpha\beta} \langle \alpha|\hat{h}_0|\beta\rangle\delta_{\alpha \alpha_1} \delta_{\beta \alpha_2} = \langle \alpha_1|\hat{h}_0|\alpha_2\rangle \label{eq:2-30f} \tag{42}
# \end{equation}
# $$

# ## Two-body operators in second quantization
# 
# Let us now derive the expression for our two-body interaction part, which also conserves the number of particles.
# We can proceed in exactly the same way as for the one-body operator. In the coordinate representation our
# two-body interaction part takes the following expression

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-31"></div>
# 
# $$
# \begin{equation}
# 	\hat{H}_I = \sum_{i < j} V(x_i,x_j) \label{eq:2-31} \tag{43}
# \end{equation}
# $$

# where the summation runs over distinct pairs. The term $V$ can be an interaction model for the nucleon-nucleon interaction
# or the interaction between two electrons. It can also include additional two-body interaction terms. 
# 
# The action of this operator on a product of 
# two single-particle functions is defined as

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-32"></div>
# 
# $$
# \begin{equation}
# 	V(x_i,x_j) \psi_{\alpha_k}(x_i) \psi_{\alpha_l}(x_j) = \sum_{\alpha_k'\alpha_l'} 
# 		\psi_{\alpha_k}'(x_i)\psi_{\alpha_l}'(x_j) 
# 		\langle \alpha_k'\alpha_l'|\hat{v}|\alpha_k\alpha_l\rangle \label{eq:2-32} \tag{44}
# \end{equation}
# $$

# We can now let $\hat{H}_I$ act on all terms in the linear combination for $|\alpha_1\alpha_2\dots\alpha_n\rangle$. Without any permutations we have

# $$
# \left( \sum_{i < j} V(x_i,x_j) \right) \psi_{\alpha_1}(x_1)\psi_{\alpha_2}(x_2)\dots \psi_{\alpha_n}(x_n) \nonumber
# $$

# $$
# = \sum_{\alpha_1'\alpha_2'} \langle \alpha_1'\alpha_2'|\hat{v}|\alpha_1\alpha_2\rangle
# 		\psi_{\alpha_1}'(x_1)\psi_{\alpha_2}'(x_2)\dots \psi_{\alpha_n}(x_n) \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# $$
# + \sum_{\alpha_1'\alpha_n'} \langle \alpha_1'\alpha_n'|\hat{v}|\alpha_1\alpha_n\rangle
# 		\psi_{\alpha_1}'(x_1)\psi_{\alpha_2}(x_2)\dots \psi_{\alpha_n}'(x_n) \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# $$
# + \sum_{\alpha_2'\alpha_n'} \langle \alpha_2'\alpha_n'|\hat{v}|\alpha_2\alpha_n\rangle
# 		\psi_{\alpha_1}(x_1)\psi_{\alpha_2}'(x_2)\dots \psi_{\alpha_n}'(x_n) \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-33"></div>
# 
# $$
# \begin{equation} 
# 	  + \dots \label{eq:2-33} \tag{45}
# \end{equation}
# $$

# where on the rhs we have a term for each distinct pairs. 
# 
# For the other terms on the rhs we obtain similar expressions  and summing over all terms we obtain

# $$
# H_I |\alpha_1\alpha_2\dots\alpha_n\rangle = \sum_{\alpha_1', \alpha_2'} \langle \alpha_1'\alpha_2'|\hat{v}|\alpha_1\alpha_2\rangle
# 		|\alpha_1'\alpha_2'\dots\alpha_n\rangle \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# $$
# + \sum_{\alpha_1', \alpha_n'} \langle \alpha_1'\alpha_n'|\hat{v}|\alpha_1\alpha_n\rangle
# 		|\alpha_1'\alpha_2\dots\alpha_n'\rangle \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# $$
# + \sum_{\alpha_2', \alpha_n'} \langle \alpha_2'\alpha_n'|\hat{v}|\alpha_2\alpha_n\rangle
# 		|\alpha_1\alpha_2'\dots\alpha_n'\rangle \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-34"></div>
# 
# $$
# \begin{equation} 
# 	 + \dots \label{eq:2-34} \tag{46}
# \end{equation}
# $$

# We introduce second quantization via the relation

# $$
# a_{\alpha_k'}^{\dagger} a_{\alpha_l'}^{\dagger} a_{\alpha_l} a_{\alpha_k} 
# 		|\alpha_1\alpha_2\dots\alpha_k\dots\alpha_l\dots\alpha_n\rangle \nonumber
# $$

# $$
# = (-1)^{k-1} (-1)^{l-2} a_{\alpha_k'}^{\dagger} a_{\alpha_l'}^{\dagger} a_{\alpha_l} a_{\alpha_k}
# 		|\alpha_k\alpha_l \underbrace{\alpha_1\alpha_2\dots\alpha_n}_{\neq \alpha_k,\alpha_l}\rangle \nonumber
# $$

# $$
# = (-1)^{k-1} (-1)^{l-2} 
# 	|\alpha_k'\alpha_l' \underbrace{\alpha_1\alpha_2\dots\alpha_n}_{\neq \alpha_k',\alpha_l'}\rangle \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-35"></div>
# 
# $$
# \begin{equation} 
# 	= |\alpha_1\alpha_2\dots\alpha_k'\dots\alpha_l'\dots\alpha_n\rangle \label{eq:2-35} \tag{47}
# \end{equation}
# $$

# Inserting this in ([46](#eq:2-34)) gives

# $$
# H_I |\alpha_1\alpha_2\dots\alpha_n\rangle
# 	= \sum_{\alpha_1', \alpha_2'} \langle \alpha_1'\alpha_2'|\hat{v}|\alpha_1\alpha_2\rangle
# 		a_{\alpha_1'}^{\dagger} a_{\alpha_2'}^{\dagger} a_{\alpha_2} a_{\alpha_1}
# 		|\alpha_1\alpha_2\dots\alpha_n\rangle \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# $$
# = \sum_{\alpha_1', \alpha_n'} \langle \alpha_1'\alpha_n'|\hat{v}|\alpha_1\alpha_n\rangle
# 		a_{\alpha_1'}^{\dagger} a_{\alpha_n'}^{\dagger} a_{\alpha_n} a_{\alpha_1}
# 		|\alpha_1\alpha_2\dots\alpha_n\rangle \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# $$
# = \sum_{\alpha_2', \alpha_n'} \langle \alpha_2'\alpha_n'|\hat{v}|\alpha_2\alpha_n\rangle
# 		a_{\alpha_2'}^{\dagger} a_{\alpha_n'}^{\dagger} a_{\alpha_n} a_{\alpha_2}
# 		|\alpha_1\alpha_2\dots\alpha_n\rangle \nonumber
# $$

# $$
# + \dots \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-36"></div>
# 
# $$
# \begin{equation} 
# 	= \sum_{\alpha, \beta, \gamma, \delta} ' \langle \alpha\beta|\hat{v}|\gamma\delta\rangle
# 		a^{\dagger}_\alpha a^{\dagger}_\beta a_\delta a_\gamma
# 		|\alpha_1\alpha_2\dots\alpha_n\rangle \label{eq:2-36} \tag{48}
# \end{equation}
# $$

# Here we let $\sum'$ indicate that the sums running over $\alpha$ and $\beta$ run over all
# single-particle states, while the summations  $\gamma$ and $\delta$ 
# run over all pairs of single-particle states. We wish to remove this restriction and since

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-37"></div>
# 
# $$
# \begin{equation}
# 	\langle \alpha\beta|\hat{v}|\gamma\delta\rangle = \langle \beta\alpha|\hat{v}|\delta\gamma\rangle \label{eq:2-37} \tag{49}
# \end{equation}
# $$

# we get

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-38a"></div>
# 
# $$
# \begin{equation}
# 	\sum_{\alpha\beta} \langle \alpha\beta|\hat{v}|\gamma\delta\rangle a^{\dagger}_\alpha a^{\dagger}_\beta a_\delta a_\gamma = 
# 		\sum_{\alpha\beta} \langle \beta\alpha|\hat{v}|\delta\gamma\rangle 
# 		a^{\dagger}_\alpha a^{\dagger}_\beta a_\delta a_\gamma \label{eq:2-38a} \tag{50} 
# \end{equation}
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-38b"></div>
# 
# $$
# \begin{equation} 
# 	= \sum_{\alpha\beta}\langle \beta\alpha|\hat{v}|\delta\gamma\rangle
# 		a^{\dagger}_\beta a^{\dagger}_\alpha a_\gamma a_\delta \label{eq:2-38b} \tag{51}
# \end{equation}
# $$

# where we  have used the anti-commutation rules.
# 
# Changing the summation indices 
# $\alpha$ and $\beta$ in ([51](#eq:2-38b)) we obtain

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-38c"></div>
# 
# $$
# \begin{equation}
# 	\sum_{\alpha\beta} \langle \alpha\beta|\hat{v}|\gamma\delta\rangle a^{\dagger}_\alpha a^{\dagger}_\beta a_\delta a_\gamma =
# 		 \sum_{\alpha\beta} \langle \alpha\beta|\hat{v}|\delta\gamma\rangle 
# 		  a^{\dagger}_\alpha a^{\dagger}_\beta  a_\gamma a_\delta \label{eq:2-38c} \tag{52}
# \end{equation}
# $$

# From this it follows that the restriction on the summation over $\gamma$ and $\delta$ can be removed if we multiply with a factor $\frac{1}{2}$, resulting in

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-39"></div>
# 
# $$
# \begin{equation}
# 	\hat{H}_I = \frac{1}{2} \sum_{\alpha\beta\gamma\delta} \langle \alpha\beta|\hat{v}|\gamma\delta\rangle
# 		a^{\dagger}_\alpha a^{\dagger}_\beta a_\delta a_\gamma \label{eq:2-39} \tag{53}
# \end{equation}
# $$

# where we sum freely over all single-particle states $\alpha$, 
# $\beta$, $\gamma$ og $\delta$.
# 
# With this expression we can now verify that the second quantization form of $\hat{H}_I$ in Eq. ([53](#eq:2-39)) 
# results in the same matrix between two anti-symmetrized two-particle states as its corresponding coordinate
# space representation. We have

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-40"></div>
# 
# $$
# \begin{equation}
# 	\langle \alpha_1 \alpha_2|\hat{H}_I|\beta_1 \beta_2\rangle =
# 		\frac{1}{2} \sum_{\alpha\beta\gamma\delta}
# 			\langle \alpha\beta|\hat{v}|\gamma\delta\rangle \langle 0|a_{\alpha_2} a_{\alpha_1} 
# 			 a^{\dagger}_\alpha a^{\dagger}_\beta a_\delta a_\gamma 
# 			 a_{\beta_1}^{\dagger} a_{\beta_2}^{\dagger}|0\rangle. \label{eq:2-40} \tag{54}
# \end{equation}
# $$

# Using the commutation relations we get

# $$
# a_{\alpha_2} a_{\alpha_1}a^{\dagger}_\alpha a^{\dagger}_\beta 
# 		a_\delta a_\gamma a_{\beta_1}^{\dagger} a_{\beta_2}^{\dagger} \nonumber
# $$

# $$
# = a_{\alpha_2} a_{\alpha_1}a^{\dagger}_\alpha a^{\dagger}_\beta 
# 		( a_\delta \delta_{\gamma \beta_1} a_{\beta_2}^{\dagger} - 
# 		a_\delta  a_{\beta_1}^{\dagger} a_\gamma a_{\beta_2}^{\dagger} ) \nonumber
# $$

# $$
# = a_{\alpha_2} a_{\alpha_1}a^{\dagger}_\alpha a^{\dagger}_\beta 
# 		(\delta_{\gamma \beta_1} \delta_{\delta \beta_2} - \delta_{\gamma \beta_1} a_{\beta_2}^{\dagger} a_\delta -
# 		a_\delta a_{\beta_1}^{\dagger}\delta_{\gamma \beta_2} +
# 		a_\delta a_{\beta_1}^{\dagger} a_{\beta_2}^{\dagger} a_\gamma ) \nonumber
# $$

# $$
# = a_{\alpha_2} a_{\alpha_1}a^{\dagger}_\alpha a^{\dagger}_\beta 
# 		(\delta_{\gamma \beta_1} \delta_{\delta \beta_2} - \delta_{\gamma \beta_1} a_{\beta_2}^{\dagger} a_\delta \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-41"></div>
# 
# $$
# \begin{equation} 
# 		 \qquad - \delta_{\delta \beta_1} \delta_{\gamma \beta_2} + \delta_{\gamma \beta_2} a_{\beta_1}^{\dagger} a_\delta
# 		+ a_\delta a_{\beta_1}^{\dagger} a_{\beta_2}^{\dagger} a_\gamma ) \label{eq:2-41} \tag{55}
# \end{equation}
# $$

# The vacuum expectation value of this product of operators becomes

# $$
# \langle 0|a_{\alpha_2} a_{\alpha_1} a^{\dagger}_\alpha a^{\dagger}_\beta a_\delta a_\gamma 
# 		a_{\beta_1}^{\dagger} a_{\beta_2}^{\dagger}|0\rangle \nonumber
# $$

# $$
# = (\delta_{\gamma \beta_1} \delta_{\delta \beta_2} -
# 		\delta_{\delta \beta_1} \delta_{\gamma \beta_2} ) 
# 		\langle 0|a_{\alpha_2} a_{\alpha_1}a^{\dagger}_\alpha a^{\dagger}_\beta|0\rangle \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-42b"></div>
# 
# $$
# \begin{equation} 
# 	= (\delta_{\gamma \beta_1} \delta_{\delta \beta_2} -\delta_{\delta \beta_1} \delta_{\gamma \beta_2} )
# 	(\delta_{\alpha \alpha_1} \delta_{\beta \alpha_2} -\delta_{\beta \alpha_1} \delta_{\alpha \alpha_2} ) \label{eq:2-42b} \tag{56}
# \end{equation}
# $$

# Insertion of 
# Eq. ([56](#eq:2-42b)) in Eq. ([54](#eq:2-40)) results in

# $$
# \langle \alpha_1\alpha_2|\hat{H}_I|\beta_1\beta_2\rangle = \frac{1}{2} \big[ 
# 		\langle \alpha_1\alpha_2|\hat{v}|\beta_1\beta_2\rangle- \langle \alpha_1\alpha_2|\hat{v}|\beta_2\beta_1\rangle \nonumber
# $$

# $$
# - \langle \alpha_2\alpha_1|\hat{v}|\beta_1\beta_2\rangle + \langle \alpha_2\alpha_1|\hat{v}|\beta_2\beta_1\rangle \big] \nonumber
# $$

# $$
# = \langle \alpha_1\alpha_2|\hat{v}|\beta_1\beta_2\rangle - \langle \alpha_1\alpha_2|\hat{v}|\beta_2\beta_1\rangle \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-43b"></div>
# 
# $$
# \begin{equation} 
# 	= \langle \alpha_1\alpha_2|\hat{v}|\beta_1\beta_2\rangle_{\mathrm{AS}}. \label{eq:2-43b} \tag{57}
# \end{equation}
# $$

# The two-body operator can also be expressed in terms of the anti-symmetrized matrix elements we discussed previously as

# $$
# \hat{H}_I = \frac{1}{2} \sum_{\alpha\beta\gamma\delta}  \langle \alpha \beta|\hat{v}|\gamma \delta\rangle
# 		a_\alpha^{\dagger} a_\beta^{\dagger} a_\delta a_\gamma \nonumber
# $$

# $$
# = \frac{1}{4} \sum_{\alpha\beta\gamma\delta} \left[ \langle \alpha \beta|\hat{v}|\gamma \delta\rangle -
# 		\langle \alpha \beta|\hat{v}|\delta\gamma \rangle \right] 
# 		a_\alpha^{\dagger} a_\beta^{\dagger} a_\delta a_\gamma \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-45"></div>
# 
# $$
# \begin{equation} 
# 	= \frac{1}{4} \sum_{\alpha\beta\gamma\delta} \langle \alpha \beta|\hat{v}|\gamma \delta\rangle_{\mathrm{AS}}
# 		a_\alpha^{\dagger} a_\beta^{\dagger} a_\delta a_\gamma \label{eq:2-45} \tag{58}
# \end{equation}
# $$

# The factors in front of the operator, either  $\frac{1}{4}$ or 
# $\frac{1}{2}$ tells whether we use antisymmetrized matrix elements or not. 
# 
# We can now express the Hamiltonian operator for a many-fermion system  in the occupation basis representation
# as

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-46b"></div>
# 
# $$
# \begin{equation}
# 	H = \sum_{\alpha, \beta} \langle \alpha|\hat{t}+\hat{u}_{\mathrm{ext}}|\beta\rangle a_\alpha^{\dagger} a_\beta + \frac{1}{4} \sum_{\alpha\beta\gamma\delta}
# 		\langle \alpha \beta|\hat{v}|\gamma \delta\rangle a_\alpha^{\dagger} a_\beta^{\dagger} a_\delta a_\gamma. \label{eq:2-46b} \tag{59}
# \end{equation}
# $$

# This is the form we will use in the rest of these lectures, assuming that we work with anti-symmetrized two-body matrix elements.

# ## Proof of Wick's theorem
# 
# Discuss also Wick's generalized theorem

# ## Interaction, Schroedinger and Heisenberg pictures

# ## Time dependent wick's theorem

# ## Gell-Man and Low's theorem and Adiabatic switching

# ## Particle-hole formalism
# 
# Second quantization is a useful and elegant formalism  for constructing many-body  states and 
# quantum mechanical operators. One can express and translate many physical processes
# into simple pictures such as Feynman diagrams. Expecation values of many-body states are also easily calculated.
# However, although the equations are seemingly easy to set up, from  a practical point of view, that is
# the solution of Schroedinger's equation, there is no particular gain.
# The many-body equation is equally hard to solve, irrespective of representation. 
# The cliche that 
# there is no free lunch brings us down to earth again.  
# Note however that a transformation to a particular
# basis, for cases where the interaction obeys specific symmetries, can ease the solution of Schroedinger's equation. 
# 
# But there is at least one important case where second quantization comes to our rescue.
# It is namely easy to introduce another reference state than the pure vacuum $|0\rangle $, where all single-particle states are active.
# With many particles present it is often useful to introduce another reference state  than the vacuum state$|0\rangle $. We will label this state $|c\rangle$ ($c$ for core) and as we will see it can reduce 
# considerably the complexity and thereby the dimensionality of the many-body problem. It allows us to sum up to infinite order specific many-body correlations.  The particle-hole representation is one of these handy representations. 
# 
# In the original particle representation these states are products of the creation operators  $a_{\alpha_i}^\dagger$ acting on the true vacuum $|0\rangle $.
# Following Eq. ([3](#eq:2-2)) we have

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-47a"></div>
# 
# $$
# \begin{equation}
#  |\alpha_1\alpha_2\dots\alpha_{n-1}\alpha_n\rangle = a_{\alpha_1}^\dagger a_{\alpha_2}^\dagger \dots
# 					a_{\alpha_{n-1}}^\dagger a_{\alpha_n}^\dagger |0\rangle  \label{eq:2-47a} \tag{60} 
# \end{equation}
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-47b"></div>
# 
# $$
# \begin{equation} 
# 	|\alpha_1\alpha_2\dots\alpha_{n-1}\alpha_n\alpha_{n+1}\rangle =
# 		a_{\alpha_1}^\dagger a_{\alpha_2}^\dagger \dots a_{\alpha_{n-1}}^\dagger a_{\alpha_n}^\dagger
# 		a_{\alpha_{n+1}}^\dagger |0\rangle  \label{eq:2-47b} \tag{61} 
# \end{equation}
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-47c"></div>
# 
# $$
# \begin{equation} 
# 	|\alpha_1\alpha_2\dots\alpha_{n-1}\rangle = a_{\alpha_1}^\dagger a_{\alpha_2}^\dagger \dots
# 		a_{\alpha_{n-1}}^\dagger |0\rangle  \label{eq:2-47c} \tag{62}
# \end{equation}
# $$

# If we use Eq. ([60](#eq:2-47a)) as our new reference state, we can simplify considerably the representation of 
# this state

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-48a"></div>
# 
# $$
# \begin{equation}
# 	|c\rangle  \equiv |\alpha_1\alpha_2\dots\alpha_{n-1}\alpha_n\rangle =
# 		a_{\alpha_1}^\dagger a_{\alpha_2}^\dagger \dots a_{\alpha_{n-1}}^\dagger a_{\alpha_n}^\dagger |0\rangle  \label{eq:2-48a} \tag{63}
# \end{equation}
# $$

# The new reference states for the $n+1$ and $n-1$ states can then be written as

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-48b"></div>
# 
# $$
# \begin{equation}
# 	|\alpha_1\alpha_2\dots\alpha_{n-1}\alpha_n\alpha_{n+1}\rangle = (-1)^n a_{\alpha_{n+1}}^\dagger |c\rangle 
# 		\equiv (-1)^n |\alpha_{n+1}\rangle_c \label{eq:2-48b} \tag{64} 
# \end{equation}
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-48c"></div>
# 
# $$
# \begin{equation} 
# 	|\alpha_1\alpha_2\dots\alpha_{n-1}\rangle = (-1)^{n-1} a_{\alpha_n} |c\rangle  
# 		\equiv (-1)^{n-1} |\alpha_{n-1}\rangle_c \label{eq:2-48c} \tag{65} 
# \end{equation}
# $$

# The first state has one additional particle with respect to the new vacuum state
# $|c\rangle $  and is normally referred to as a one-particle state or one particle added to the 
# many-body reference state. 
# The second state has one particle less than the reference vacuum state  $|c\rangle $ and is referred to as
# a one-hole state. 
# When dealing with a new reference state it is often convenient to introduce 
# new creation and annihilation operators since we have 
# from Eq. ([65](#eq:2-48c))

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-49"></div>
# 
# $$
# \begin{equation}
# 	a_\alpha |c\rangle  \neq 0 \label{eq:2-49} \tag{66}
# \end{equation}
# $$

# since  $\alpha$ is contained  in $|c\rangle $, while for the true vacuum we have 
# $a_\alpha |0\rangle  = 0$ for all $\alpha$.
# 
# The new reference state leads to the definition of new creation and annihilation operators
# which satisfy the following relations

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-50a"></div>
# 
# $$
# \begin{equation}
# 	b_\alpha |c\rangle  = 0 \label{eq:2-50a} \tag{67} 
# \end{equation}
# $$

# $$
# \{b_\alpha^\dagger , b_\beta^\dagger \} = \{b_\alpha , b_\beta \} = 0 \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-50c"></div>
# 
# $$
# \begin{equation} 
# 	\{b_\alpha^\dagger , b_\beta \} = \delta_{\alpha \beta} \label{eq:2-50c} \tag{68}
# \end{equation}
# $$

# We assume also that the new reference state is properly normalized

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-51"></div>
# 
# $$
# \begin{equation}
# 	\langle c | c \rangle = 1 \label{eq:2-51} \tag{69}
# \end{equation}
# $$

# The physical interpretation of these new operators is that of so-called quasiparticle states.
# This means that a state defined by the addition of one extra particle to a reference state $|c\rangle $ may not necesseraly be interpreted as one particle coupled to a core.
# We define now new creation operators that act on a state $\alpha$ creating a new quasiparticle state

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-52"></div>
# 
# $$
# \begin{equation}
# 	b_\alpha^\dagger|c\rangle  = \Bigg\{ \begin{array}{ll}
# 		a_\alpha^\dagger |c\rangle  = |\alpha\rangle, & \alpha > F \\
# 		\\
# 		a_\alpha |c\rangle  = |\alpha^{-1}\rangle, & \alpha \leq F
# 	\end{array} \label{eq:2-52} \tag{70}
# \end{equation}
# $$

# where $F$ is the Fermi level representing the last  occupied single-particle orbit 
# of the new reference state $|c\rangle $. 
# 
# The annihilation is the hermitian conjugate of the creation operator

# $$
# b_\alpha = (b_\alpha^\dagger)^\dagger,
# $$

# resulting in

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-54"></div>
# 
# $$
# \begin{equation}
# 	b_\alpha^\dagger = \Bigg\{ \begin{array}{ll}
# 		a_\alpha^\dagger & \alpha > F \\
# 		\\
# 		a_\alpha & \alpha \leq F
# 	\end{array} \qquad 
# 	b_\alpha = \Bigg\{ \begin{array}{ll}
# 		a_\alpha & \alpha > F \\
# 		\\
# 		 a_\alpha^\dagger & \alpha \leq F
# 	\end{array} \label{eq:2-54} \tag{71}
# \end{equation}
# $$

# With the new creation and annihilation operator  we can now construct 
# many-body quasiparticle states, with one-particle-one-hole states, two-particle-two-hole
# states etc in the same fashion as we previously constructed many-particle states. 
# We can write a general particle-hole state as

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-56"></div>
# 
# $$
# \begin{equation}
# 	|\beta_1\beta_2\dots \beta_{n_p} \gamma_1^{-1} \gamma_2^{-1} \dots \gamma_{n_h}^{-1}\rangle \equiv
# 		\underbrace{b_{\beta_1}^\dagger b_{\beta_2}^\dagger \dots b_{\beta_{n_p}}^\dagger}_{>F}
# 		\underbrace{b_{\gamma_1}^\dagger b_{\gamma_2}^\dagger \dots b_{\gamma_{n_h}}^\dagger}_{\leq F} |c\rangle \label{eq:2-56} \tag{72}
# \end{equation}
# $$

# We can now rewrite our one-body and two-body operators in terms of the new creation and annihilation operators.
# The number operator becomes

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-57b"></div>
# 
# $$
# \begin{equation}
# 	\hat{N} = \sum_\alpha a_\alpha^\dagger a_\alpha= 
# \sum_{\alpha > F} b_\alpha^\dagger b_\alpha + n_c - \sum_{\alpha \leq F} b_\alpha^\dagger b_\alpha \label{eq:2-57b} \tag{73}
# \end{equation}
# $$

# where $n_c$ is the number of particle in the new vacuum state $|c\rangle $.  
# The action of $\hat{N}$ on a many-body state results in

# <!-- Equation labels as ordinary links -->
# <div id="2-59"></div>
# 
# $$
# \begin{equation}
# 	N |\beta_1\beta_2\dots \beta_{n_p} \gamma_1^{-1} \gamma_2^{-1} \dots \gamma_{n_h}^{-1}\rangle = (n_p + n_c - n_h) |\beta_1\beta_2\dots \beta_{n_p} \gamma_1^{-1} \gamma_2^{-1} \dots \gamma_{n_h}^{-1}\rangle \label{2-59} \tag{74}
# \end{equation}
# $$

# Here  $n=n_p +n_c - n_h$ is the total number of particles in the quasi-particle state of 
# Eq. ([72](#eq:2-56)). Note that  $\hat{N}$ counts the total number of particles  present

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-60"></div>
# 
# $$
# \begin{equation}
# 	N_{qp} = \sum_\alpha b_\alpha^\dagger b_\alpha, \label{eq:2-60} \tag{75}
# \end{equation}
# $$

# gives us the number of quasi-particles as can be seen by computing

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-61"></div>
# 
# $$
# \begin{equation}
# 	N_{qp}= |\beta_1\beta_2\dots \beta_{n_p} \gamma_1^{-1} \gamma_2^{-1} \dots \gamma_{n_h}^{-1}\rangle
# 		= (n_p + n_h)|\beta_1\beta_2\dots \beta_{n_p} \gamma_1^{-1} \gamma_2^{-1} \dots \gamma_{n_h}^{-1}\rangle \label{eq:2-61} \tag{76}
# \end{equation}
# $$

# where $n_{qp} = n_p + n_h$ is the total number of quasi-particles.
# 
# We express the one-body operator $\hat{H}_0$ in terms of the quasi-particle creation and annihilation operators, resulting in

# $$
# \hat{H}_0 = \sum_{\alpha\beta > F} \langle \alpha|\hat{h}_0|\beta\rangle  b_\alpha^\dagger b_\beta +
# 		\sum_{\alpha > F, \beta \leq F } \left[\langle \alpha|\hat{h}_0|\beta\rangle b_\alpha^\dagger b_\beta^\dagger + \langle \beta|\hat{h}_0|\alpha\rangle b_\beta  b_\alpha \right] \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-63b"></div>
# 
# $$
# \begin{equation} 
# 	+ \sum_{\alpha \leq F} \langle \alpha|\hat{h}_0|\alpha\rangle - \sum_{\alpha\beta \leq F} \langle \beta|\hat{h}_0|\alpha\rangle b_\alpha^\dagger b_\beta \label{eq:2-63b} \tag{77}
# \end{equation}
# $$

# The first term  gives contribution only for particle states, while the last one
# contributes only for holestates. The second term can create or destroy a set of
# quasi-particles and 
# the third term is the contribution  from the vacuum state $|c\rangle$.
# 
# Before we continue with the expressions for the two-body operator, we introduce a nomenclature we will use for the rest of this
# text. It is inspired by the notation used in quantum chemistry.
# We reserve the labels $i,j,k,\dots$ for hole states and $a,b,c,\dots$ for states above $F$, viz. particle states.
# This means also that we will skip the constraint $\leq F$ or $> F$ in the summation symbols. 
# Our operator $\hat{H}_0$  reads now

# $$
# \hat{H}_0 = \sum_{ab} \langle a|\hat{h}|b\rangle b_a^\dagger b_b +
# 		\sum_{ai} \left[
# 		\langle a|\hat{h}|i\rangle b_a^\dagger b_i^\dagger + 
# 		\langle i|\hat{h}|a\rangle b_i  b_a \right] \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-63c"></div>
# 
# $$
# \begin{equation} 
# 	+ \sum_{i} \langle i|\hat{h}|i\rangle - 
# 		\sum_{ij} \langle j|\hat{h}|i\rangle
# 		b_i^\dagger b_j \label{eq:2-63c} \tag{78}
# \end{equation}
# $$

# The two-particle operator in the particle-hole formalism  is more complicated since we have
# to translate four indices $\alpha\beta\gamma\delta$ to the possible combinations of particle and hole
# states.  When performing the commutator algebra we can regroup the operator in five different terms

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-65"></div>
# 
# $$
# \begin{equation}
# 	\hat{H}_I = \hat{H}_I^{(a)} + \hat{H}_I^{(b)} + \hat{H}_I^{(c)} + \hat{H}_I^{(d)} + \hat{H}_I^{(e)} \label{eq:2-65} \tag{79}
# \end{equation}
# $$

# Using anti-symmetrized  matrix elements, 
# bthe term  $\hat{H}_I^{(a)}$ is

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-66"></div>
# 
# $$
# \begin{equation}
# 	\hat{H}_I^{(a)} = \frac{1}{4}
# 	\sum_{abcd} \langle ab|\hat{V}|cd\rangle 
# 		b_a^\dagger b_b^\dagger b_d b_c \label{eq:2-66} \tag{80}
# \end{equation}
# $$

# The next term $\hat{H}_I^{(b)}$  reads

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-67b"></div>
# 
# $$
# \begin{equation}
# 	 \hat{H}_I^{(b)} = \frac{1}{4} \sum_{abci}\left(\langle ab|\hat{V}|ci\rangle b_a^\dagger b_b^\dagger b_i^\dagger b_c +\langle ai|\hat{V}|cb\rangle b_a^\dagger b_i b_b b_c\right) \label{eq:2-67b} \tag{81}
# \end{equation}
# $$

# This term conserves the number of quasiparticles but creates or removes a 
# three-particle-one-hole  state. 
# For $\hat{H}_I^{(c)}$  we have

# $$
# \hat{H}_I^{(c)} = \frac{1}{4}
# 		\sum_{abij}\left(\langle ab|\hat{V}|ij\rangle b_a^\dagger b_b^\dagger b_j^\dagger b_i^\dagger +
# 		\langle ij|\hat{V}|ab\rangle b_a  b_b b_j b_i \right)+  \nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-68c"></div>
# 
# $$
# \begin{equation} 
# 		\frac{1}{2}\sum_{abij}\langle ai|\hat{V}|bj\rangle b_a^\dagger b_j^\dagger b_b b_i + 
# 		\frac{1}{2}\sum_{abi}\langle ai|\hat{V}|bi\rangle b_a^\dagger b_b. \label{eq:2-68c} \tag{82}
# \end{equation}
# $$

# The first line stands for the creation of a two-particle-two-hole state, while the second line represents
# the creation to two one-particle-one-hole pairs
# while the last term represents a contribution to the particle single-particle energy
# from the hole states, that is an interaction between the particle states and the hole states
# within the new vacuum  state.
# The fourth term reads

# $$
# \hat{H}_I^{(d)} = \frac{1}{4} 
# 	 	\sum_{aijk}\left(\langle ai|\hat{V}|jk\rangle b_a^\dagger b_k^\dagger b_j^\dagger b_i+
# \langle ji|\hat{V}|ak\rangle b_k^\dagger b_j b_i b_a\right)+\nonumber
# $$

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-69d"></div>
# 
# $$
# \begin{equation} 
# \frac{1}{4}\sum_{aij}\left(\langle ai|\hat{V}|ji\rangle b_a^\dagger b_j^\dagger+
# \langle ji|\hat{V}|ai\rangle - \langle ji|\hat{V}|ia\rangle b_j b_a \right). \label{eq:2-69d} \tag{83} 
# \end{equation}
# $$

# The terms in the first line  stand for the creation of a particle-hole state 
# interacting with hole states, we will label this as a two-hole-one-particle contribution. 
# The remaining terms are a particle-hole state interacting with the holes in the vacuum state. 
# Finally we have

# <!-- Equation labels as ordinary links -->
# <div id="eq:2-70d"></div>
# 
# $$
# \begin{equation}
# 	\hat{H}_I^{(e)} = \frac{1}{4}
# 		 \sum_{ijkl}
# 		 \langle kl|\hat{V}|ij\rangle b_i^\dagger b_j^\dagger b_l b_k+
# 	        \frac{1}{2}\sum_{ijk}\langle ij|\hat{V}|kj\rangle b_k^\dagger b_i
# 	        +\frac{1}{2}\sum_{ij}\langle ij|\hat{V}|ij\rangle \label{eq:2-70d} \tag{84}
# \end{equation}
# $$

# The first terms represents the 
# interaction between two holes while the second stands for the interaction between a hole and the remaining holes in the vacuum state.
# It represents a contribution to single-hole energy  to first order.
# The last term collects all contributions to the energy of the ground state of a closed-shell system arising
# from hole-hole correlations.

# ## Summarizing and defining a normal-ordered Hamiltonian

# $$
# \Phi_{AS}(\alpha_1, \dots, \alpha_A; x_1, \dots x_A)=
#             \frac{1}{\sqrt{A}} \sum_{\hat{P}} (-1)^P \hat{P} \prod_{i=1}^A \psi_{\alpha_i}(x_i),
# $$

# which is equivalent with $|\alpha_1 \dots \alpha_A\rangle= a_{\alpha_1}^{\dagger} \dots a_{\alpha_A}^{\dagger} |0\rangle$. We have also

# $$
# a_p^\dagger|0\rangle = |p\rangle, \quad a_p |q\rangle = \delta_{pq}|0\rangle
# $$

# $$
# \delta_{pq} = \left\{a_p, a_q^\dagger \right\},
# $$

# and

# $$
# 0 = \left\{a_p^\dagger, a_q \right\} = \left\{a_p, a_q \right\} = \left\{a_p^\dagger, a_q^\dagger \right\}
# $$

# $$
# |\Phi_0\rangle = |\alpha_1 \dots \alpha_A\rangle, \quad \alpha_1, \dots, \alpha_A \leq \alpha_F
# $$

# $$
# \left\{a_p^\dagger, a_q \right\}= \delta_{pq}, p, q \leq \alpha_F
# $$

# $$
# \left\{a_p, a_q^\dagger \right\} = \delta_{pq}, p, q > \alpha_F
# $$

# with         $i,j,\ldots \leq \alpha_F, \quad a,b,\ldots > \alpha_F, \quad p,q, \ldots - \textrm{any}$

# $$
# a_i|\Phi_0\rangle = |\Phi_i\rangle, \hspace{0.5cm} a_a^\dagger|\Phi_0\rangle = |\Phi^a\rangle
# $$

# and

# $$
# a_i^\dagger|\Phi_0\rangle = 0 \hspace{0.5cm}  a_a|\Phi_0\rangle = 0
# $$

# The one-body operator is defined as

# $$
# \hat{F} = \sum_{pq} \langle p|\hat{f}|q\rangle a_p^\dagger a_q
# $$

# while the two-body opreator is defined as

# $$
# \hat{V} = \frac{1}{4} \sum_{pqrs} \langle pq|\hat{v}|rs\rangle_{AS} a_p^\dagger a_q^\dagger a_s a_r
# $$

# where we have defined the antisymmetric matrix elements

# $$
# \langle pq|\hat{v}|rs\rangle_{AS} = \langle pq|\hat{v}|rs\rangle - \langle pq|\hat{v}|sr\rangle.
# $$

# We can also define a three-body operator

# $$
# \hat{V}_3 = \frac{1}{36} \sum_{pqrstu} \langle pqr|\hat{v}_3|stu\rangle_{AS} 
#                 a_p^\dagger a_q^\dagger a_r^\dagger a_u a_t a_s
# $$

# with the antisymmetrized matrix element

# <!-- Equation labels as ordinary links -->
# <div id="_auto2"></div>
# 
# $$
# \begin{equation}
#             \langle pqr|\hat{v}_3|stu\rangle_{AS} = \langle pqr|\hat{v}_3|stu\rangle + \langle pqr|\hat{v}_3|tus\rangle + \langle pqr|\hat{v}_3|ust\rangle- \langle pqr|\hat{v}_3|sut\rangle - \langle pqr|\hat{v}_3|tsu\rangle - \langle pqr|\hat{v}_3|uts\rangle.
# \label{_auto2} \tag{85}
# \end{equation}
# $$

# ## Operators in second quantization
# 
# In the build-up of an FCI code that is meant to tackle large dimensionalities
# is the action of the Hamiltonian $\hat{H}$ on a
# Slater determinant represented in second quantization as

# $$
# |\alpha_1\dots \alpha_n\rangle = a_{\alpha_1}^{\dagger} a_{\alpha_2}^{\dagger} \dots a_{\alpha_n}^{\dagger} |0\rangle.
# $$

# The time consuming part stems from the action of the Hamiltonian
# on the above determinant,

# $$
# \left(\sum_{\alpha\beta} \langle \alpha|t+u|\beta\rangle a_\alpha^{\dagger} a_\beta + \frac{1}{4} \sum_{\alpha\beta\gamma\delta}
#                 \langle \alpha \beta|\hat{v}|\gamma \delta\rangle a_\alpha^{\dagger} a_\beta^{\dagger} a_\delta a_\gamma\right)a_{\alpha_1}^{\dagger} a_{\alpha_2}^{\dagger} \dots a_{\alpha_n}^{\dagger} |0\rangle.
# $$

# A practically useful way to implement this action is to encode a Slater determinant as a bit pattern.
# 
# Assume that we have at our disposal $n$ different single-particle orbits
# $\alpha_0,\alpha_2,\dots,\alpha_{n-1}$ and that we can distribute  among these orbits $N\le n$ particles.
# 
# A Slater  determinant can then be coded as an integer of $n$ bits. As an example, if we have $n=16$ single-particle states
# $\alpha_0,\alpha_1,\dots,\alpha_{15}$ and $N=4$ fermions occupying the states $\alpha_3$, $\alpha_6$, $\alpha_{10}$ and $\alpha_{13}$
# we could write this Slater determinant as

# $$
# \Phi_{\Lambda} = a_{\alpha_3}^{\dagger} a_{\alpha_6}^{\dagger} a_{\alpha_{10}}^{\dagger} a_{\alpha_{13}}^{\dagger} |0\rangle.
# $$

# The unoccupied single-particle states have bit value $0$ while the occupied ones are represented by bit state $1$. 
# In the binary notation we would write this   16 bits long integer as

# $$
# \begin{array}{cccccccccccccccc}
# {\alpha_0}&{\alpha_1}&{\alpha_2}&{\alpha_3}&{\alpha_4}&{\alpha_5}&{\alpha_6}&{\alpha_7} & {\alpha_8} &{\alpha_9} & {\alpha_{10}} &{\alpha_{11}} &{\alpha_{12}} &{\alpha_{13}} &{\alpha_{14}} & {\alpha_{15}} \\
# {0} & {0} &{0} &{1} &{0} &{0} &{1} &{0} &{0} &{0} &{1} &{0} &{0} &{1} &{0} & {0} \\
# \end{array}
# $$

# which translates into the decimal number

# $$
# 2^3+2^6+2^{10}+2^{13}=9288.
# $$

# We can thus encode a Slater determinant as a bit pattern.
# 
# With $N$ particles that can be distributed over $n$ single-particle states, the total number of Slater determinats (and defining thereby the dimensionality of the system) is

# $$
# \mathrm{dim}(\mathcal{H}) = \left(\begin{array}{c} n \\N\end{array}\right).
# $$

# The total number of bit patterns is $2^n$. 
# 
# We assume again that we have at our disposal $n$ different single-particle orbits
# $\alpha_0,\alpha_2,\dots,\alpha_{n-1}$ and that we can distribute  among these orbits $N\le n$ particles.
# The ordering among these states is important as it defines the order of the creation operators.
# We will write the determinant

# $$
# \Phi_{\Lambda} = a_{\alpha_3}^{\dagger} a_{\alpha_6}^{\dagger} a_{\alpha_{10}}^{\dagger} a_{\alpha_{13}}^{\dagger} |0\rangle,
# $$

# in a more compact way as

# $$
# \Phi_{3,6,10,13} = |0001001000100100\rangle.
# $$

# The action of a creation operator is thus

# $$
# a^{\dagger}_{\alpha_4}\Phi_{3,6,10,13} = a^{\dagger}_{\alpha_4}|0001001000100100\rangle=a^{\dagger}_{\alpha_4}a_{\alpha_3}^{\dagger} a_{\alpha_6}^{\dagger} a_{\alpha_{10}}^{\dagger} a_{\alpha_{13}}^{\dagger} |0\rangle,
# $$

# which becomes

# $$
# -a_{\alpha_3}^{\dagger} a^{\dagger}_{\alpha_4} a_{\alpha_6}^{\dagger} a_{\alpha_{10}}^{\dagger} a_{\alpha_{13}}^{\dagger} |0\rangle=-|0001101000100100\rangle.
# $$

# Similarly

# $$
# a^{\dagger}_{\alpha_6}\Phi_{3,6,10,13} = a^{\dagger}_{\alpha_6}|0001001000100100\rangle=a^{\dagger}_{\alpha_6}a_{\alpha_3}^{\dagger} a_{\alpha_6}^{\dagger} a_{\alpha_{10}}^{\dagger} a_{\alpha_{13}}^{\dagger} |0\rangle,
# $$

# which becomes

# $$
# -a^{\dagger}_{\alpha_4} (a_{\alpha_6}^{\dagger})^ 2 a_{\alpha_{10}}^{\dagger} a_{\alpha_{13}}^{\dagger} |0\rangle=0!
# $$

# This gives a simple recipe:  
# * If one of the bits $b_j$ is $1$ and we act with a creation operator on this bit, we return a null vector
# 
# * If $b_j=0$, we set it to $1$ and return a sign factor $(-1)^l$, where $l$ is the number of bits set before bit $j$.
# 
# Consider the action of $a^{\dagger}_{\alpha_2}$ on various slater determinants:

# $$
# \begin{array}{ccc}
# a^{\dagger}_{\alpha_2}\Phi_{00111}& = a^{\dagger}_{\alpha_2}|00111\rangle&=0\times |00111\rangle\\
# a^{\dagger}_{\alpha_2}\Phi_{01011}& = a^{\dagger}_{\alpha_2}|01011\rangle&=(-1)\times |01111\rangle\\
# a^{\dagger}_{\alpha_2}\Phi_{01101}& = a^{\dagger}_{\alpha_2}|01101\rangle&=0\times |01101\rangle\\
# a^{\dagger}_{\alpha_2}\Phi_{01110}& = a^{\dagger}_{\alpha_2}|01110\rangle&=0\times |01110\rangle\\
# a^{\dagger}_{\alpha_2}\Phi_{10011}& = a^{\dagger}_{\alpha_2}|10011\rangle&=(-1)\times |10111\rangle\\
# a^{\dagger}_{\alpha_2}\Phi_{10101}& = a^{\dagger}_{\alpha_2}|10101\rangle&=0\times |10101\rangle\\
# a^{\dagger}_{\alpha_2}\Phi_{10110}& = a^{\dagger}_{\alpha_2}|10110\rangle&=0\times |10110\rangle\\
# a^{\dagger}_{\alpha_2}\Phi_{11001}& = a^{\dagger}_{\alpha_2}|11001\rangle&=(+1)\times |11101\rangle\\
# a^{\dagger}_{\alpha_2}\Phi_{11010}& = a^{\dagger}_{\alpha_2}|11010\rangle&=(+1)\times |11110\rangle\\
# \end{array}
# $$

# What is the simplest way to obtain the phase when we act with one annihilation(creation) operator
# on the given Slater determinant representation?
# 
# We have an SD representation

# $$
# \Phi_{\Lambda} = a_{\alpha_0}^{\dagger} a_{\alpha_3}^{\dagger} a_{\alpha_6}^{\dagger} a_{\alpha_{10}}^{\dagger} a_{\alpha_{13}}^{\dagger} |0\rangle,
# $$

# in a more compact way as

# $$
# \Phi_{0,3,6,10,13} = |1001001000100100\rangle.
# $$

# The action of

# $$
# a^{\dagger}_{\alpha_4}a_{\alpha_0}\Phi_{0,3,6,10,13} = a^{\dagger}_{\alpha_4}|0001001000100100\rangle=a^{\dagger}_{\alpha_4}a_{\alpha_3}^{\dagger} a_{\alpha_6}^{\dagger} a_{\alpha_{10}}^{\dagger} a_{\alpha_{13}}^{\dagger} |0\rangle,
# $$

# which becomes

# $$
# -a_{\alpha_3}^{\dagger} a^{\dagger}_{\alpha_4} a_{\alpha_6}^{\dagger} a_{\alpha_{10}}^{\dagger} a_{\alpha_{13}}^{\dagger} |0\rangle=-|0001101000100100\rangle.
# $$

# The action

# $$
# a_{\alpha_0}\Phi_{0,3,6,10,13} = |0001001000100100\rangle,
# $$

# can be obtained by subtracting the logical sum (AND operation) of $\Phi_{0,3,6,10,13}$ and 
# a word which represents only $\alpha_0$, that is

# $$
# |1000000000000000\rangle,
# $$

# from $\Phi_{0,3,6,10,13}= |1001001000100100\rangle$.
# 
# This operation gives $|0001001000100100\rangle$. 
# 
# Similarly, we can form $a^{\dagger}_{\alpha_4}a_{\alpha_0}\Phi_{0,3,6,10,13}$, say, by adding 
# $|0000100000000000\rangle$ to $a_{\alpha_0}\Phi_{0,3,6,10,13}$, first checking that their logical sum
# is zero in order to make sure that orbital $\alpha_4$ is not already occupied. 
# 
# It is trickier however to get the phase $(-1)^l$. 
# One possibility is as follows
# * Let $S_1$ be a word that represents the $1-$bit to be removed and all others set to zero.
# 
# In the previous example $S_1=|1000000000000000\rangle$
# * Define $S_2$ as the similar word that represents the bit to be added, that is in our case
# 
# $S_2=|0000100000000000\rangle$.
# * Compute then $S=S_1-S_2$, which here becomes

# $$
# S=|0111000000000000\rangle
# $$

# * Perform then the logical AND operation of $S$ with the word containing

# $$
# \Phi_{0,3,6,10,13} = |1001001000100100\rangle,
# $$

# which results in $|0001000000000000\rangle$. Counting the number of $1-$bits gives the phase.  Here you need however an algorithm for bitcounting. Several efficient ones available.

# ## Schr\"odinger picture
# 
# The time-dependent Schr\"odinger equation (or equation of motion) reads

# $$
# \imath \hbar\frac{\partial }{\partial t}|\Psi_S(t)\rangle = \hat{H}\Psi_S(t)\rangle,
# $$

# where the subscript $S$ stands for Schr\"odinger here.
# A formal solution is given by

# $$
# |\Psi_S(t)\rangle = \exp{(-\imath\hat{H}(t-t_0)/\hbar)}|\Psi_S(t_0)\rangle.
# $$

# The Hamiltonian $\hat{H}$ is hermitian and the exponent represents a unitary 
# operator with an operation carried ut on the wave function at a time $t_0$.

# ## Interaction picture
# Our Hamiltonian is normally written out as the sum of an unperturbed part $\hat{H}_0$ and an interaction part $\hat{H}_I$, that is

# $$
# \hat{H}=\hat{H}_0+\hat{H}_I.
# $$

# In general we have $[\hat{H}_0,\hat{H}_I]\ne 0$ since $[\hat{T},\hat{V}]\ne 0$.
# We wish now to define a unitary transformation in terms of $\hat{H}_0$ by defining

# $$
# |\Psi_I(t)\rangle = \exp{(\imath\hat{H}_0t/\hbar)}|\Psi_S(t)\rangle,
# $$

# which is again a unitary transformation carried out now at the time $t$ on the 
# wave function in the Schr\"odinger picture. 
# 
# We can easily find the equation of motion by taking the time derivative

# $$
# \imath \hbar\frac{\partial }{\partial t}|\Psi_I(t)\rangle = -\hat{H}_0\exp{(\imath\hat{H}_0t/\hbar)}\Psi_S(t)\rangle+\exp{(\imath\hat{H}_0t/\hbar)}
# \imath \hbar\frac{\partial }{\partial t}\Psi_S(t)\rangle.
# $$

# Using the definition of the Schr\"odinger equation, we can rewrite the last equation as

# $$
# \imath \hbar\frac{\partial }{\partial t}|\Psi_I(t)\rangle = \exp{(\imath\hat{H}_0t/\hbar)}\left[-\hat{H}_0+\hat{H}_0+\hat{H}_I\right]\exp{(-\imath\hat{H}_0t/\hbar)}\Psi_I(t)\rangle,
# $$

# which gives us

# $$
# \imath \hbar\frac{\partial }{\partial t}|\Psi_I(t)\rangle = \hat{H}_I(t)\Psi_I(t)\rangle,
# $$

# with

# $$
# \hat{H}_I(t)=
# \exp{(\imath\hat{H}_0t/\hbar)}\hat{H}_I\exp{(-\imath\hat{H}_0t/\hbar)}.
# $$

# The order of the operators is important since $\hat{H}_0$ and $\hat{H}_I$ do generally not commute.
# The expectation value of
# an arbitrary operator in the interaction picture can now be written as

# $$
# \langle \Psi'_S(t)|\hat{O}_S|\Psi_S(t)\rangle = 
# \langle \Psi'_I(t) |\exp{(\imath\hat{H}_0t/\hbar)}\hat{O}_I
# \exp{(-\imath\hat{H}_0t/\hbar)}|\Psi_I(t)\rangle,
# $$

# and using the definition

# $$
# \hat{O}_I(t)=
# \exp{(\imath\hat{H}_0t/\hbar)}\hat{O}_I\exp{(-\imath\hat{H}_0t/\hbar)},
# $$

# we obtain

# $$
# \langle \Psi'_S(t)|\hat{O}_S|\Psi_S(t)\rangle = 
# \langle \Psi'_I(t) |\hat{O}_I(t)|\Psi_I(t)\rangle,
# $$

# stating that a unitary transformation does not change expectation values!
# 
# If we take the time derivative of the operator in the interaction picture we arrive at the following equation of motion

# $$
# \imath \hbar\frac{\partial }{\partial t}\hat{O}_I(t) = \exp{(\imath\hat{H}_0t/\hbar)}\left[\hat{O}_S\hat{H}_0-\hat{H}_0\hat{O}_S\right]\exp{(-\imath\hat{H}_0t/\hbar)}=\left[\hat{O}_I(t),\hat{H}_0\right].
# $$

# Here we have used the time-independence of the Schr\"odinger equation
# together with the observation that any function of an operator commutes with the operator itself. 
# 
# In order to solve the equation of motion equation in the interaction picture, we define a unitary operator
# time-development operator $\hat{U}(t,t')$. Later we will derive its
# connection with the linked-diagram theorem, which yields a
# linked expression for the actual operator. 
# The action of the operator on the wave function is

# $$
# |\Psi_I(t) \rangle = \hat{U}(t,t_0)|\Psi_I(t_0)\rangle,
# $$

# with the obvious value $\hat{U}(t_0,t_0)=1$.
# 
# The time-development operator $U$ has the
# properties that

# $$
# \hat{U}^{\dagger}(t,t')\hat{U}(t,t')=\hat{U}(t,t')\hat{U}^{\dagger}(t,t')=1,
# $$

# which implies that $U$ is unitary

# $$
# \hat{U}^{\dagger}(t,t')=\hat{U}^{-1}(t,t').
# $$

# Further,

# $$
# \hat{U}(t,t')\hat{U}(t't'')=\hat{U}(t,t'')
# $$

# and

# $$
# \hat{U}(t,t')\hat{U}(t',t)=1,
# $$

# which leads to

# $$
# \hat{U}(t,t')=\hat{U}^{\dagger}(t',t).
# $$

# Using our definition of Schr\"odinger's equation in the interaction picture, we can then construct the operator $\hat{U}$. We have defined

# $$
# |\Psi_I(t)\rangle = \exp{(\imath\hat{H}_0t/\hbar)}|\Psi_S(t)\rangle,
# $$

# which can be rewritten as

# $$
# |\Psi_I(t)\rangle = \exp{(\imath\hat{H}_0t/\hbar)}\exp{(-\imath\hat{H}(t-t_0)/\hbar)}|\Psi_S(t_0)\rangle,
# $$

# or

# $$
# |\Psi_I(t)\rangle = \exp{(\imath\hat{H}_0t/\hbar)}\exp{(-\imath\hat{H}(t-t_0)/\hbar)}\exp{(-\imath\hat{H}_0t_0/\hbar)}|\Psi_I(t_0)\rangle.
# $$

# From the last expression we can define

# $$
# \hat{U}(t,t_0)=\exp{(\imath\hat{H}_0t/\hbar)}\exp{(-\imath\hat{H}(t-t_0)/\hbar)}\exp{(-\imath\hat{H}_0t_0/\hbar)}.
# $$

# It is then easy to convince oneself that the properties defined above are satisfied by the definition of $\hat{U}$. 
# 
# We derive the equation of motion for $\hat{U}$ using the above definition.
# This results in

# $$
# \imath \hbar\frac{\partial }{\partial t}\hat{U}(t,t_0) = \hat{H}_I(t)\hat{U}(t,t_0),
# $$

# which we integrate from $t_0$ to a time $t$ resulting in

# $$
# \hat{U}(t,t_0)-\hat{U}(t_0,t_0)=\hat{U}(t,t_0)-1=-\frac{\imath}{\hbar}\int_{t_0}^t dt' \hat{H}_I(t')\hat{U}(t',t_0),
# $$

# which can be rewritten as

# $$
# \hat{U}(t,t_0)=1-\frac{\imath}{\hbar}\int_{t_0}^t dt' \hat{H}_I(t')\hat{U}(t',t_0).
# $$

# We can solve this equation iteratively keeping in mind the time-ordering of the of the operators

# $$
# \hat{U}(t,t_0)=1-\frac{\imath}{\hbar}\int_{t_0}^t dt' \hat{H}_I(t')+\left(\frac{-\imath}{\hbar}\right)^2\int_{t_0}^t dt'\int_{t_0}^{t'} dt'' \hat{H}_I(t')\hat{H}_I(t'')+\dots
# $$

# The third term can be written as

# $$
# \int_{t_0}^t dt'\int_{t_0}^{t'} dt'' \hat{H}_I(t')\hat{H}_I(t'')=
# \frac{1}{2}\int_{t_0}^t dt'\int_{t_0}^{t'} dt'' \hat{H}_I(t')\hat{H}_I(t'')
# +\frac{1}{2}\int_{t_0}^t dt''\int_{t''}^{t} dt' \hat{H}_I(t')\hat{H}_I(t'').
# $$

# We obtain this expression by changing the integration order in the second term
# via a change of the integration variables $t'$ and $t''$  in

# $$
# \frac{1}{2}\int_{t_0}^t dt'\int_{t_0}^{t'} dt'' \hat{H}_I(t')\hat{H}_I(t'').
# $$

# We can rewrite the terms which contain the double integral as

# $$
# \int_{t_0}^t dt'\int_{t_0}^{t'} dt'' \hat{H}_I(t')\hat{H}_I(t'')=
# $$

# $$
# \frac{1}{2}\int_{t_0}^t dt'\int_{t_0}^{t'} dt''\left[\hat{H}_I(t')\hat{H}_I(t'')\Theta(t'-t'')
# +\hat{H}_I(t')\hat{H}_I(t'')\Theta(t''-t')\right],
# $$

# with $\Theta(t''-t')$ being the standard Heavyside or step function. The step function allows us to give a specific time-ordering to the above expression.
# 
# With the $\Theta$-function we can rewrite the last expression as

# $$
# \int_{t_0}^t dt'\int_{t_0}^{t'} dt'' \hat{H}_I(t')\hat{H}_I(t'')=
# \frac{1}{2}\int_{t_0}^t dt'\int_{t_0}^{t'} dt''\hat{T}\left[\hat{H}_I(t')\hat{H}_I(t'')\right],
# $$

# where $\Hat{T}$ is the so-called time-ordering operator. 
# 
# With this definition, we can rewrite the expression for $\hat{U}$ as

# $$
# \hat{U}(t,t_0)=\sum_{n=0}^{\infty}\left(\frac{-\imath}{\hbar}\right)^n\frac{1}{n1}
# \int_{t_0}^t dt_1\dots \int_{t_0}^t dt_N \hat{T}\left[\hat{H}_I(t_1)\dots\hat{H}_I(t_n)\right]=\hat{T}\exp{\left[\frac{-\imath}{\hbar}
# \int_{t_0}^t dt' \hat{H}_I(t')\right]}.
# $$

# The above time-evolution operator in the interaction picture will be used
# to derive various contributions to many-body perturbation theory.

# ## Heisenberg picture
# 
# We wish now to define a unitary transformation in terms of $\hat{H}$ by defining

# $$
# |\Psi_H(t)\rangle = \exp{(\imath\hat{H}t/\hbar)}|\Psi_S(t)\rangle,
# $$

# which is again a unitary transformation carried out now at the time $t$ on the 
# wave function in the Schr\"odinger picture. If we combine this equation with 
# Schr\"odinger's equation we obtain the following equation of motion

# $$
# \imath \hbar\frac{\partial }{\partial t}|\Psi_H(t)\rangle = 0,
# $$

# meaning that $|\Psi_H(t)\rangle$ is time independent. An operator in this picture is defined as

# $$
# \hat{O}_H(t)=
# \exp{(\imath\hat{H}t/\hbar)}\hat{O}_S\exp{(-\imath\hat{H}t/\hbar)}.
# $$

# The time dependence is then in the operator itself, and this yields in turn the
# following equation of motion

# $$
# \imath \hbar\frac{\partial }{\partial t}\hat{O}_H(t) = \exp{(\imath\hat{H}t/\hbar)}\left[\hat{O}_H\hat{H}-\hat{H}\hat{O}_H\right]\exp{(-\imath\hat{H}t/\hbar)}=\left[\hat{O}_H(t),\hat{H}\right].
# $$

# We note that an operator in the Heisenberg picture can be related to the corresponding
# operator in the interaction picture as

# $$
# \hat{O}_H(t)=
# \exp{(\imath\hat{H}t/\hbar)}\hat{O}_S\exp{(-\imath\hat{H}t/\hbar)}=
# 
# \exp{(\imath\hat{H}_It/\hbar)}\exp{(-\imath\hat{H}_0t/\hbar)}\hat{O}_I\exp{(\imath\hat{H}_0t/\hbar)}\exp{(-\imath\hat{H}_It/\hbar)}.
# $$

# With our definition of the time evolution operator we see that

# $$
# \hat{O}_H(t)=\hat{U}(0,t)\hat{O}_I\hat{U}(t,0),
# $$

# which in turn implies that $\hat{O}_S=\hat{O}_I(0)=\hat{O}_H(0)$, all operators are equal at $t=0$. The wave function in the Heisenberg formalism is 
# related to the other pictures as

# $$
# |\Psi_H\rangle=|\Psi_S(0)\rangle=|\Psi_I(0)\rangle,
# $$

# since the wave function in the Heisenberg picture is time independent. 
# We can relate this wave function to that a given time $t$ via the time evolution operator as

# $$
# |\Psi_H\rangle=\hat{U}(0,t)|\Psi_I(t)\rangle.
# $$

# ## Adiabatic hypothesis
# 
# We assume that the interaction term is switched on gradually. Our wave function at time $t=-\infty$ and $t=\infty$ is supposed to represent a non-interacting system
# given by the solution to the unperturbed part of our Hamiltonian $\hat{H}_0$.
# We assume the ground state is given by $|\Phi_0\rangle$, which could be a Slater determinant.
# 
# We define our Hamiltonian as

# $$
# \hat{H}=\hat{H}_0+\exp{(-\varepsilon t/\hbar)}\hat{H}_I,
# $$

# where $\varepsilon$ is a small number. The way we write the Hamiltonian 
# and its interaction term is meant to simulate the switching of the interaction.
# 
# The time evolution of the wave function in the interaction picture is then

# $$
# |\Psi_I(t) \rangle = \hat{U}_{\varepsilon}(t,t_0)|\Psi_I(t_0)\rangle,
# $$

# with

# $$
# \hat{U}_{\varepsilon}(t,t_0)=\sum_{n=0}^{\infty}\left(\frac{-\imath}{\hbar}\right)^n\frac{1}{n!}
# \int_{t_0}^t dt_1\dots \int_{t_0}^t dt_N \exp{(-\varepsilon(t_1+\dots+t_n)/\hbar)}\hat{T}\left[\hat{H}_I(t_1)\dots\hat{H}_I(t_n)\right]
# $$

# In the limit $t_0\rightarrow -\infty$, the solution ot Schr\"odinger's equation is
# $|\Phi_0\rangle$, and the eigenenergies are given by

# $$
# \hat{H}_0|\Phi_0\rangle=W_0|\Phi_0\rangle,
# $$

# meaning that

# $$
# |\Psi_S(t_0)\rangle = \exp{(-\imath W_0t_0/\hbar)}|\Phi_0\rangle,
# $$

# with the corresponding interaction picture wave function given by

# $$
# |\Psi_I(t_0)\rangle = \exp{(\imath \hat{H}_0t_0/\hbar)}|\Psi_S(t_0)\rangle=|\Phi_0\rangle.
# $$

# The solution becomes time independent in the limit $t_0\rightarrow -\infty$.
# The same conclusion can be reached by looking at

# $$
# \imath \hbar\frac{\partial }{\partial t}|\Psi_I(t)\rangle =
# \exp{(-\varepsilon |t|/\hbar)}\hat{H}_I|\Psi_I(t)\rangle
# $$

# and taking the limit $t\rightarrow \pm\infty$.
# We can rewrite the equation for the wave function at a time $t=0$ as

# $$
# |\Psi_I(0) \rangle = \hat{U}_{\varepsilon}(0,-\infty)|\Phi_0\rangle.
# $$

# ## Wigner-Jordan transformation and second quantization

# ## Baker-Campbell-Hausdorf
# Discuss also Suzuki-Trotter as an approximation to BCH

# ## Exercise 1: Relation between basis functions
# 
# This exercise serves to convince you about the relation between
# two different single-particle bases, where one could be our new Hartree-Fock basis and the other a harmonic oscillator basis.
# 
# Consider a Slater determinant built up of single-particle orbitals $\psi_{\lambda}$, 
# with $\lambda = 1,2,\dots,A$. The unitary transformation

# $$
# \psi_a  = \sum_{\lambda} C_{a\lambda}\phi_{\lambda},
# $$

# brings us into the new basis.  
# The new basis has quantum numbers $a=1,2,\dots,A$.
# Show that the new basis is orthonormal.
# Show that the new Slater determinant constructed from the new single-particle wave functions can be
# written as the determinant based on the previous basis and the determinant of the matrix $C$.
# Show that the old and the new Slater determinants are equal up to a complex constant with absolute value unity.
# (Hint, $C$ is a unitary matrix). 
# 
# Starting with the second quantization representation of the Slater determinant

# $$
# \Phi_{0}=\prod_{i=1}^{n}a_{\alpha_{i}}^{\dagger}|0\rangle,
# $$

# use Wick's theorem to compute the normalization integral
# $\langle\Phi_{0}|\Phi_{0}\rangle$.

# ## Exercise 2: Matrix elements
# 
# Calculate the matrix elements

# $$
# \langle \alpha_{1}\alpha_{2}|\hat{F}|\alpha_{1}\alpha_{2}\rangle
# $$

# and

# $$
# \langle \alpha_{1}\alpha_{2}|\hat{G}|\alpha_{1}\alpha_{2}\rangle
# $$

# with

# $$
# |\alpha_{1}\alpha_{2}\rangle=a_{\alpha_{1}}^{\dagger}a_{\alpha_{2}}^{\dagger}|0\rangle ,
# $$

# $$
# \hat{F}=\sum_{\alpha\beta}\langle \alpha|\hat{f}|\beta\rangle
# a_{\alpha}^{\dagger}a_{\beta}  ,
# $$

# $$
# \langle \alpha|\hat{f}|\beta\rangle=\int \psi_{\alpha}^{*}(x)f(x)\psi_{\beta}(x)dx ,
# $$

# $$
# \hat{G} = \frac{1}{2}\sum_{\alpha\beta\gamma\delta}
# \langle \alpha\beta |\hat{g}|\gamma\delta\rangle
# a_{\alpha}^{\dagger}a_{\beta}^{\dagger}a_{\delta}a_{\gamma} ,
# $$

# and

# $$
# \langle \alpha\beta |\hat{g}|\gamma\delta\rangle=
# \int\int \psi_{\alpha}^{*}(x_{1})\psi_{\beta}^{*}(x_{2})g(x_{1},
# x_{2})\psi_{\gamma}(x_{1})\psi_{\delta}(x_{2})dx_{1}dx_{2}
# $$

# Compare these results with those from exercise 3c).

# ## Exercise 3: Normal-ordered one-body operator
# 
# Show that the onebody part of the Hamiltonian

# $$
# \hat{H}_0 = \sum_{pq} \langle p|\hat{h}_0|q\rangle a^{\dagger}_p a_q,
# $$

# can be written, using standard annihilation and creation operators, in normal-ordered form as

# $$
# \hat{H}_0 = \sum_{pq} \langle p|\hat{h}_0|q\rangle \left\{a^\dagger_p a_q\right\} +
#              \sum_i \langle i|\hat{h}_0|i\rangle.
# $$

# Explain the meaning of the various symbols. Which reference 
# vacuum has been used?

# ## Exercise 4: Normal-ordered two-body operator
# 
# Show that the twobody part of the Hamiltonian

# $$
# \hat{H}_I = \frac{1}{4} \sum_{pqrs} \langle pq|\hat{v}|rs\rangle a^\dagger_p a^\dagger_q a_s  a_r,
# $$

# can be written, using standard annihilation and creation operators, in normal-ordered form as

# $$
# \hat{H}_I =\frac{1}{4} \sum_{pqrs} \langle pq|\hat{v}|rs\rangle \left\{a^\dagger_p a^\dagger_q a_s  a_r\right\}
#             + \sum_{pqi} \langle pi|\hat{v}|qi\rangle \left\{a^\dagger_p a_q\right\} 
#             + \frac{1}{2} \sum_{ij}\langle ij|\hat{v}|ij\rangle.
# $$

# Explain again the meaning of the various symbols.
# 
# This exercise is optional: Derive the normal-ordered form of the threebody part of the Hamiltonian.

# $$
# \hat{H}_3 = \frac{1}{36} \sum_{\substack{pqr \\ stu}}
#                  \langle pqr|\hat{v}_3|stu\rangle a^\dagger_p a^\dagger_q a^\dagger_r a_u a_t a_s,
# $$

# and specify the contributions to the twobody, onebody and the scalar part.

# ## Exercise 5: Matrix elements using the Slater-Condon rule
# 
# The aim of this exercise is to set up specific matrix elements that will turn useful when we start our discussions of the nuclear shell model. In particular you will notice, depending on the character of the operator, that many matrix elements will actually be zero.
# 
# Consider three $N$-particle  Slater determinants  $|\Phi_0$, $|\Phi_i^a\rangle$ and $|\Phi_{ij}^{ab}\rangle$, where the notation means that 
# Slater determinant $|\Phi_i^a\rangle$ differs from $|\Phi_0\rangle$ by one single-particle state, that is a single-particle
# state $\psi_i$ is replaced by a single-particle state $\psi_a$. 
# It is often interpreted as a so-called one-particle-one-hole excitation.
# Similarly, the Slater determinant $|\Phi_{ij}^{ab}\rangle$
# differs by two single-particle states from $|\Phi_0\rangle$ and is normally thought of as a two-particle-two-hole excitation.
# We assume also that $|\Phi_0\rangle$ represents our new vacuum reference state
# and the labels $ijk\dots$ represent single-particle states below the Fermi level and $abc\dots$ represent states above the Fermi level, so-called particle states.
# We define thereafter a general onebody normal-ordered (with respect to the new vacuum state) operator 
# as

# $$
# \hat{F}_N=\sum_{pq}\langle p |f |\beta\rangle \left\{a_{p}^{\dagger}a_{q}\right\}  ,
# $$

# with

# $$
# \langle p |f| q\rangle=\int \psi_{p}^{*}(x)f(x)\psi_{q}(x)dx ,
# $$

# and a general normal-ordered two-body operator

# $$
# \hat{G}_N = \frac{1}{4}\sum_{pqrs}
# \langle pq |g| rs\rangle_{AS} \left\{a_{p}^{\dagger}a_{q}^{\dagger}a_{s}a_{r}\right\} ,
# $$

# with for example the direct matrix element given as

# $$
# \langle pq |g| rs\rangle=
# \int\int \psi_{p}^{*}(x_{1})\psi_{q}^{*}(x_{2})g(x_{1}, x_{2})\psi_{r}(x_{1})\psi_{s}(x_{2})dx_{1}dx_{2}
# $$

# with $g$ being invariant under the interchange of the coordinates of two particles.
# The single-particle states $\psi_i$ are not necessarily eigenstates of $\hat{f}$.  The curly brackets mean that the operators are normal-ordered with respect to the new vacuum reference state. 
# 
# How would you write the above Slater determinants in a second quantization formalism, utilizing the fact that we have defined a new reference state? 
# 
# Use thereafter Wick's theorem to find the expectation values of

# $$
# \langle \Phi_0 \vert\hat{F}_N\vert\Phi_0\rangle,
# $$

# and

# $$
# \langle \Phi_0\hat{G}_N|\Phi_0\rangle.
# $$

# Find thereafter

# $$
# \langle \Phi_0 |\hat{F}_N|\Phi_i^a\rangle,
# $$

# and

# $$
# \langle \Phi_0|\hat{G}_N|\Phi_i^a\rangle,
# $$

# Finally, find

# $$
# \langle \Phi_0 |\hat{F}_N|\Phi_{ij}^{ab}\rangle,
# $$

# and

# $$
# \langle \Phi_0|\hat{G}_N|\Phi_{ij}^{ab}\rangle.
# $$

# What happens with the two-body operator if we have a transition probability  of the type

# $$
# \langle \Phi_0|\hat{G}_N|\Phi_{ijk}^{abc}\rangle,
# $$

# where the Slater determinant to the right of the operator differs by more than two single-particle states?

# ## Exercise 6: Program to set up Slater determinants
# 
# Write a program which sets up all possible Slater determinants given $N=4$ eletrons which can occupy
# the atomic single-particle states defined by the $1s$, $2s2p$ and $3s3p3d$ shells. How many single-particle
# states $n$ are there in total?  Include the spin degrees of freedom as well.

# ## Exercise 7: Using sympy to compute matrix elements
# 
# Compute the matrix element

# $$
# \langle\alpha_{1}\alpha_{2}\alpha_{3}|\hat{G}|\alpha_{1}'\alpha_{2}'\alpha_{3}'\rangle,
# $$

# using Wick's theorem and express the two-body operator
# $G$ in the occupation number (second quantization) 
# representation.

# ## Exercise 8: Using sympy to compute matrix elements
# 
# The last exercise can be solved using the symbolic Python package called *SymPy*. SymPy is a Python 
# package for general purpose symbolic algebra. There is a physics module with several interesting submodules.
# Among these, the submodule called *secondquant*, contains several  functionalities that allow us to test
# our algebraic manipulations using Wick's theorem and operators for second quantization.

# In[1]:


from sympy import *
from sympy.physics.secondquant import *

i, j = symbols('i,j', below_fermi=True)
a, b = symbols('a,b', above_fermi=True)
p, q = symbols('p,q')
print simplify(wicks(Fd(i)*F(a)*Fd(p)*F(q)*Fd(b)*F(j), keep_only_fully_contracted=True))


# The code defines single-particle states above and below the Fermi level, in addition to the genereal symbols
# $pq$ which can refer to any type of state below or above the Fermi level. Wick's theorem is implemented between 
# the creation and annihilation operators *Fd* and *F*, respectively. Using the simplify option, one can lump together several Kronecker-$\delta$ functions.

# ## Exercise 9: Using sympy to compute matrix elements
# 
# We can expand the above Python code by defining one-body and two-body operators using  the following SymPy code

# In[2]:


# This code sets up a two-body Hamiltonian for fermions
from sympy import symbols, latex, WildFunction, collect, Rational
from sympy.physics.secondquant import F, Fd, wicks, AntiSymmetricTensor, substitute_dummies, NO

# setup hamiltonian
p,q,r,s = symbols('p q r s',dummy=True)
f = AntiSymmetricTensor('f',(p,),(q,))
pr = NO((Fd(p)*F(q)))
v = AntiSymmetricTensor('v',(p,q),(r,s))
pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))
Hamiltonian=f*pr + Rational(1)/Rational(4)*v*pqsr
print "Hamiltonian defined as:", latex(Hamiltonian)


# Here we have used the *AntiSymmetricTensor* functionality, together with normal-ordering defined by the *NO* function. 
# Using the *latex* option, this program produces the following output

# $$
# f^{p}_{q} \left\{a^\dagger_{p} a_{q}\right\} - \frac{1}{4} v^{qp}_{sr} \left\{a^\dagger_{p} a^\dagger_{q} a_{r} a_{s}\right\}
# $$

# ## Exercise 10: Using sympy to compute matrix elements
# 
# We can now use this code to compute the matrix elements between two two-body Slater determinants using Wick's theorem.

# In[3]:


from sympy import symbols, latex, WildFunction, collect, Rational, simplify
from sympy.physics.secondquant import F, Fd, wicks, AntiSymmetricTensor, substitute_dummies, NO, evaluate_deltas
# setup hamiltonian
p,q,r,s = symbols('p q r s',dummy=True)
f = AntiSymmetricTensor('f',(p,),(q,))
pr = NO((Fd(p)*F(q)))
v = AntiSymmetricTensor('v',(p,q),(r,s))
pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))
Hamiltonian=f*pr + Rational(1)/Rational(4)*v*pqsr
c,d = symbols('c, d',above_fermi=True)
a,b = symbols('a, b',above_fermi=True)

expression = wicks(F(b)*F(a)*Hamiltonian*Fd(c)*Fd(d),keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
expression = evaluate_deltas(expression)
expression = simplify(expression)
print "Hamiltonian defined as:", latex(expression)


# The result is as expected,

# $$
# \delta_{a c} f^{b}_{d} - \delta_{a d} f^{b}_{c} - \delta_{b c} f^{a}_{d} + \delta_{b d} f^{a}_{c} + v^{ab}_{cd}.
# $$

# ## Exercise 11: Using sympy to compute matrix elements
# 
# We can continue along these lines and define a normal-ordered Hamiltonian with respect to a given reference state.
# In our first step we just define the Hamiltonian

# In[4]:


from sympy import symbols, latex, WildFunction, collect, Rational, simplify
from sympy.physics.secondquant import F, Fd, wicks, AntiSymmetricTensor, substitute_dummies, NO, evaluate_deltas
# setup hamiltonian
p,q,r,s = symbols('p q r s',dummy=True)
f = AntiSymmetricTensor('f',(p,),(q,))
pr = Fd(p)*F(q)
v = AntiSymmetricTensor('v',(p,q),(r,s))
pqsr = Fd(p)*Fd(q)*F(s)*F(r)
#define the Hamiltonian
Hamiltonian = f*pr + Rational(1)/Rational(4)*v*pqsr
#define indices for states above and below the Fermi level
index_rule = {
     'below':  'kl',
     'above':  'cd',
     'general': 'pqrs'
     }
Hnormal = substitute_dummies(Hamiltonian,new_indices=True, pretty_indices=index_rule)
print "Hamiltonian defined as:", latex(Hnormal)


# which results in

# $$
# f^{q}_{p} a^\dagger_{q} a_{p} + \frac{1}{4} v^{sr}_{qp} a^\dagger_{s} a^\dagger_{r} a_{p} a_{q}
# $$

# ## Exercise 12: Using sympy to compute matrix elements
# 
# In our next step we define the reference energy $E_0$ and redefine the Hamiltonian by subtracting the reference energy and collecting the coefficients for all normal-ordered products (given by the *NO* function).

# In[5]:


from sympy import symbols, latex, WildFunction, collect, Rational, simplify
from sympy.physics.secondquant import F, Fd, wicks, AntiSymmetricTensor, substitute_dummies, NO, evaluate_deltas
# setup hamiltonian
p,q,r,s = symbols('p q r s',dummy=True)
f = AntiSymmetricTensor('f',(p,),(q,))
pr = Fd(p)*F(q)
v = AntiSymmetricTensor('v',(p,q),(r,s))
pqsr = Fd(p)*Fd(q)*F(s)*F(r)
#define the Hamiltonian
Hamiltonian=f*pr + Rational(1)/Rational(4)*v*pqsr
#define indices for states above and below the Fermi level
index_rule = {
     'below':  'kl',
     'above':  'cd',
     'general': 'pqrs'
     }
Hnormal = substitute_dummies(Hamiltonian,new_indices=True, pretty_indices=index_rule)
E0 = wicks(Hnormal,keep_only_fully_contracted=True)
Hnormal = Hnormal-E0
w = WildFunction('w')
Hnormal = collect(Hnormal, NO(w))
Hnormal = evaluate_deltas(Hnormal)
print latex(Hnormal)


# which gives us

# $$
# - f^{i}_{i} + f^{q}_{p} a^\dagger_{q} a_{p} - \frac{1}{4} v^{ii}_{ii} - \frac{1}{4} v^{ii}_{ii} + \frac{1}{4} v^{sr}_{qp} a^\dagger_{r} a^\dagger_{s} a_{q} a_{p},
# $$

# again as expected, with the reference energy to be subtracted.

# ## Exercise 13: Using sympy to compute matrix elements
# 
# We can now go back to exercise 7 and define the Hamiltonian and the second-quantized representation of a  three-body Slater determinant.

# In[6]:


from sympy import symbols, latex, WildFunction, collect, Rational, simplify
from sympy.physics.secondquant import F, Fd, wicks, AntiSymmetricTensor, substitute_dummies, NO, evaluate_deltas
# setup hamiltonian
p,q,r,s = symbols('p q r s',dummy=True)

v = AntiSymmetricTensor('v',(p,q),(r,s))
pqsr = NO(Fd(p)*Fd(q)*F(s)*F(r))
Hamiltonian=Rational(1)/Rational(4)*v*pqsr
a,b,c,d,e,f = symbols('a,b, c, d, e, f',above_fermi=True)

expression = wicks(F(c)*F(b)*F(a)*Hamiltonian*Fd(d)*Fd(e)*Fd(f),keep_only_fully_contracted=True, simplify_kronecker_deltas=True)
expression = evaluate_deltas(expression)
expression = simplify(expression)
print latex(expression)


# resulting in nine terms (as expected),

# $$
# - \delta_{a d} v^{cb}_{ef} - \delta_{a e} v^{cb}_{fd} + \delta_{a f} v^{cb}_{ed} - \delta_{b d} v^{ac}_{ef} - \delta_{b e} v^{ac}_{fd} + \delta_{b f} v^{ac}_{ed} + \delta_{c d} v^{ab}_{ef} + \delta_{c e} v^{ab}_{fd} - \delta_{c f} v^{ab}_{ed}
# $$

# ## Exercise 14: Diagrammatic representation of Hartree-Fock equations
# 
# What is the diagrammatic representation of the HF equation?

# $$
# -\langle\alpha_{k}|u^{HF}|\alpha_{i}\rangle+\sum_{j=1}^{n}
# \left[\langle\alpha_{k}\alpha_{j}|\hat{v}|\alpha_{i}\alpha_{j}\rangle-
# \langle\alpha_{k}\alpha_{j}|v|\alpha_{j}\alpha_{i}\rangle\right]=0
# $$

# (Represent $(-u^{HF})$ by the symbol $---$X .)

# ## Exercise 15: Derivation of Hartree-Fock equations
# 
# Consider the ground state $|\Phi\rangle$ 
# of a bound many-particle system of fermions. Assume that we remove one particle
# from the single-particle state $\lambda$ and that our system ends in a new state
# $|\Phi_{n}\rangle$. 
# Define the energy needed to remove this particle as

# $$
# E_{\lambda}=\sum_{n}\vert\langle\Phi_{n}|a_{\lambda}|\Phi\rangle\vert^{2}(E_{0}-E_{n}),
# $$

# where $E_{0}$ and $E_{n}$  are the ground state energies of the states
# $|\Phi\rangle$  and  $|\Phi_{n}\rangle$, respectively.
#  * Show that

# $$
# E_{\lambda}=\langle\Phi|a_{\lambda}^{\dagger}\left[
# a_{\lambda},H \right]|\Phi\rangle,
# $$

# where $H$ is the Hamiltonian of this system.
#  * If we assume that $\Phi$ is the  Hartree-Fock result, find the 
# 
# relation between $E_{\lambda}$ and the single-particle energy
# $\varepsilon_{\lambda}$
# for states $\lambda \leq F$ and $\lambda >F$, with

# $$
# \varepsilon_{\lambda}=\langle\lambda|\hat{t}+\hat{u}|\lambda\rangle,
# $$

# and

# $$
# \langle\lambda|\hat{u}|\lambda\rangle=\sum_{\beta \leq F}
# \langle\lambda\beta|\hat{v}|\lambda\beta\rangle.
# $$

# We have assumed an antisymmetrized matrix element here.
# Discuss the result.
# 
# The Hamiltonian operator is defined as

# $$
# H=\sum_{\alpha\beta}\langle\alpha|\hat{t}|\beta\rangle a_{\alpha}^{\dagger}a_{\beta}+
# \frac{1}{2}\sum_{\alpha\beta\gamma\delta}\langle\alpha\beta|\hat{v}|\gamma\delta\rangle a_{\alpha}^{\dagger}a_{\beta}^{\dagger}a_{\delta}a_{\gamma}.
# $$
