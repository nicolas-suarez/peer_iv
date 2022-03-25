For one of my research projects I was working on peer effects for college students, and I wanted to use the Generalized 2SLS procedure to estimate peer effects described in [Bramoullé, Y., Djebbari, H., & Fortin, B. (2009). Identification of peer effects through social networks. Journal of econometrics, 150(1), 41-55.](https://www.sciencedirect.com/science/article/pii/S0304407609000335). The problem is that this project required a lot of data work and cleaning, so I worked on it on Stata mostly, but there is not a lot of support for peer effects regressions on `Stata`, so I had to adapt the code of Bramoullé et al. (2009) on my own. I based my code on their article, and on the code shared on 
[Habiba Djebbari's website](https://rpubs.com/Nicolas_Gallo/549370).

Here I share my code, some examples about how to use it, and a comparison of my code and the original `R` code used by the authors.

To wrote this, I took advantage of Stata's new integration with Python and Jupyter Notebooks (more information on how to use this feature can be found [here](https://www.stata.com/new-in-stata/jupyter-notebooks/)). I also used the [rpy2 package](https://rpy2.github.io/) to run `R` commands inside a Jupyter notebook.

## Setting up Stata in Python


```python
#setting up Stata in Python
import stata_setup
stata_setup.config("C:/Program Files/Stata17", "mp")
```

    
      ___  ____  ____  ____  ____ ®
     /__    /   ____/   /   ____/      17.0
    ___/   /   /___/   /   /___/       MP—Parallel Edition
    
     Statistics and Data Science       Copyright 1985-2021 StataCorp LLC
                                       StataCorp
                                       4905 Lakeway Drive
                                       College Station, Texas 77845 USA
                                       800-STATA-PC        https://www.stata.com
                                       979-696-4600        stata@stata.com
    
    Stata license: Unlimited-user 4-core network, expiring 21 Jul 2022
    Serial number: xxxxxxxxxxxx
      Licensed to: Nicolas Suarez
                   Stanford University
    
    Notes:
          1. Unicode is supported; see help unicode_advice.
          2. More than 2 billion observations are allowed; see help obs_advice.
          3. Maximum number of variables is set to 5,000; see help set_maxvar.
    
    Running C:\Program Files\Stata17/profile.do ...
    

# Peer IV function

To implement linear in means regressions, I wrote the `peer_iv` function. This function takes as inputs, as any other regression function, a dependent variable and independent variables, plus an adjacency matrix `G` to perform the calculations. This matrix has to be a `Mata` matrix object (`Stata` matrices have size limitations). The option `row` allows us to row-normalize the adjacency matrix (so the sum of each row is 1, and we can interpret the product of a variable and the matrix as weighted means), and the `fixed` option adds group level fixed effects.

The matrix generates an standard Stata regression output, containing coefficients, standard errors, p-values and all the relevant information, and it stores eclass results. This means that the output could be stored with custom commands like `outreg2`, `estout` or `esttab` that allow the user to build customizable output tables.


```python
%%stata

capture program drop peer_iv
program define peer_iv, eclass
version 17
syntax varlist, ADJacency(name) [ROW FIXED OLS] 
/* implements the generalized 2SLS model of Bramoulle et al (2009), without fixed effects.
The model includes a constant, then the endogenous effect, effects of the independent variables, and then the exogenous effects.
For more details see https://rpubs.com/Nicolas_Gallo/549370

INPUTS:
varlist = exogenous variables
dep= dependent variable
adjacency=name of the adjacency matrix in Mata
row=optional, row normalizes the adjacency matrix
fixed=optional, estimates model with cohort level fixed effects
ols=optional, OLS results. Don't use together with FIXED

OUTPUT:
displays the coefficients and standard errors in a table. Stores eclass results.
*/

*separating dependent from independent variables
gettoken dep varlist: varlist
preserve
quietly{
*checking if there are missing values in our data
reg `dep' `varlist'
*recovering the indexes of non-missing observations
gen muestra=e(sample)
ereturn clear 

*moving data as matrices
mata X=st_data(.,"`varlist'")
mata y=st_data(.,"`dep'")
mata muestra=st_data(.,"muestra")

*dropping missing values from data matrices
mata X=select(X,muestra)
mata y=select(y,muestra)

*dropping missing values from G matrix (eliminating the rows and columns with missing values, so the matrix are comformable)
mata G1=select(`adjacency',muestra)
mata G1=select(G1,muestra')

*row normalizing G if needed
if "`row'"!="" mata G1=G1:/editvalue(rowsum(G1),0,1)

*generating identity matrix
mata Id=I(rows(G1))

*OLS results
if "`ols'"!="" {
	mata X_1 =  J(rows(X),1,1), G1*y, X, G1*X 
	mata theta= invsym(quadcross(X_1, X_1))*quadcross(X_1, y)
	mata e= y - X_1*theta
	mata V = (quadsum(e:^2)/(rows(X_1)-cols(X_1)))*invsym(quadcross(X_1, X_1))
}
else {
	*putting matrices together
	*with fixed effects
	if "`fixed'"!="" {
		mata S=( (Id-G1)*X, (Id-G1)*G1*X, (Id-G1)*G1*G1*X )
		mata X_1= ( (Id-G1)*G1*y, (Id-G1)*X, (Id-G1)*G1*X )				
	}
	else{
		mata S=( J(rows(X),1,1), X, G1*X, G1*G1*X )
		mata X_1= ( J(rows(X),1,1), G1*y, X, G1*X )
	}
	mata P= S*invsym(quadcross(S,S))*S'

	*first 2sls
	if "`fixed'"!="" mata theta_1= invsym(X_1'*P*X_1)*X_1'*P*(Id-G1)*y
    else mata theta_1= invsym(X_1'*P*X_1)*X_1'*P*y
	
	*building instrument
	if "`fixed'"!="" {
		mata Z = G1*luinv(Id-theta_1[1]*G1)*(Id-G1)*(X*theta_1[2::(1+cols(X))] +  G1*X*theta_1[(2+cols(X))::(1+2*cols(X))] ), (Id-G1)*X, (Id-G1)*G1*X	
	}
	else{
		mata Z = J(rows(X),1,1), G1*luinv(Id-theta_1[2]*G1)*( theta_1[1]*J(rows(X),1,1) + X*theta_1[3::(2+cols(X))] +  G1*X*theta_1[(3+cols(X))::(2+2*cols(X))] ), X, G1*X
	}
	*
	
    *final 2sls
    if "`fixed'"!="" mata theta = luinv(quadcross(Z,X_1))*quadcross(Z,(Id-G1)*y)
    else mata theta = luinv(quadcross(Z,X_1))*quadcross(Z,y)

	*resids
	if "`fixed'"!="" {
		mata e= (Id-G1)*y - luinv(Id-theta[1]*G1)*((Id-G1)*X*theta[2::(1+cols(X))] + (Id-G1)*G1*X*theta[(2+cols(X))::(1+2*cols(X))] )
	}
	else{
		mata e= y - luinv(Id-theta[2]*G1)*( theta[1]*J(rows(X),1,1) + X*theta[3::(2+cols(X))] +  G1*X*theta[(3+cols(X))::(2+2*cols(X))] )
	}


	*variance
	mata V = luinv(quadcross(Z,X_1))*(Z')*diag(e:^2)*Z*luinv(quadcross(X_1,Z))
}

*sending results to Stata
mata st_matrix("b",theta')
mata st_matrix("V",V)

*row and col names for matrices
local exog_peer //list for names of exogenous effects
foreach var in `varlist'{
	local exog_peer `exog_peer' `var'_p
}
if "`fixed'"!="" {
	local varnames `dep'_p `varlist' `exog_peer'
}
else{
	local varnames _cons `dep'_p `varlist' `exog_peer'
}


*adding col and rownames
matrix colnames b= `varnames'
matrix colnames V = `varnames'
matrix rownames V = `varnames'
}
*storing eclass results
ereturn post b V, depname(`dep') esample(muestra)
mata st_numscalar("e(N)", rows(G1))
mata st_numscalar("e(df_r)", rows(X_1)-cols(X_1))
eret local cmd peer_iv
ereturn display

restore		
end

```

## Example

Here we will run a little example to see how the command works, and how you can generate an adjacency matrix in Stata. We are going to use the `auto` dataset with 1978 automobile data, and we are going to create a random adjacency matrix `G` with elements that are drawn from a uniform between 0 an 1, but we will force the elements in the main diagonal to be 0. We are also going to row normalize the adjacency matrix.


```python
%%stata

sysuse auto, clear

*generating the matrix
mata G= runiform(`c(N)', `c(N)') 
forval i=1/`c(N)'{
	mata G[strtoreal(st_local("i")),strtoreal(st_local("i"))]=0
}

*running the regression
peer_iv price trunk turn, row adj(G)

eret list
```

    
    . 
    . sysuse auto, clear
    (1978 automobile data)
    
    . 
    . *generating the matrix
    . mata G= runiform(`c(N)', `c(N)') 
    
    . forval i=1/`c(N)'{
      2.         mata G[strtoreal(st_local("i")),strtoreal(st_local("i"))]=0
      3. }
    
    . 
    . *running the regression
    . peer_iv price trunk turn, row adj(G)
    ------------------------------------------------------------------------------
           price | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
    -------------+----------------------------------------------------------------
           _cons |   41706.23   48765.18     0.86   0.395    -55603.18    139015.6
         price_p |  -1.527064   11.64602    -0.13   0.896    -24.76633     21.7122
           trunk |   141.5504   83.60055     1.69   0.095    -25.27191    308.3727
            turn |    98.8758   103.6205     0.95   0.343    -107.8956    305.6472
         trunk_p |     439.69   968.8704     0.45   0.651    -1493.661    2373.041
          turn_p |  -959.4334   2872.551    -0.33   0.739    -6691.519    4772.652
    ------------------------------------------------------------------------------
    
    . 
    . eret list
    
    scalars:
                      e(N) =  74
                   e(df_r) =  68
    
    macros:
                    e(cmd) : "peer_iv"
             e(properties) : "b V"
                 e(depvar) : "price"
    
    matrices:
                      e(b) :  1 x 6
                      e(V) :  6 x 6
    
    . 
    

We can see the output of the regression, including a constant, the coefficient `price_p` is the coefficient of the endogenous effect, while `trunk_p` and `turn_p` are the exogenous effects.

After the regression command we ran the `eret list` command, and we can see all the elements that are stored after running the command.

## Storing Mata matrices

At least for my particular application, computing the adjacency matrix was very slow, so it is not something that I would do each time before I want to run peer effects regressions. To avoid this, we can use the `Mata` functions `matsave` and `matuse` to store a matrix as a `.mmat` object, and then load it into Stata.

# Checking if the code works
Here, to see if my code works properly, I will run the `R` code provided by Bramoullé, Djebbari and Fortin (it can be found [here](https://rpubs.com/Nicolas_Gallo/549370)) to generate data, then estimate peer effects with and without fixed effects, export their data to Stata, and then see if my function obtains the same coefficients.

The code provided by the authors is meant to be used to run Monte Carlo simulations, so I made some modifications to keep only the relevant parts. Also, for both cases, the authors defined different data generating processes, so we will have 2 vectors of dependent variables, but all of them are generated with the vector of white noise.


```python
#Python package to use magic R commands
%load_ext rpy2.ipython
```    

## Generating network data in R

```r
%%R

library(knitr)
library(igraph)
library(truncnorm)
set.seed(1)

alpha=0.7683
beta=0.4666
gamma=0.0834
delta=0.1507
e_var=0.1
#Generating a graph with 100 vertices and a probability of link of 0.04 with the "random.renyi.game()" function
g<-erdos.renyi.game(100,0.04)
#Generating the associated weighted adjacency matrix
G<-get.adjacency(g)
G<-as.matrix(G)

#Drawing a vector x of characteristics
x_sim<-matrix(rbinom(n = nrow(G),size = 1,prob = 0.9458 ),nrow(G),1)
for(i in 1:nrow(x_sim)){
  if(x_sim[i,] != 0){
    x_sim[i,]<-rtruncnorm(n = 1,a = 0,b = 1000,mean = 1,sd = 3) 
  }
}
#a vector filled with 1,  size m x 1 (used when there is an intercept in the model, i.e. when fixed_effects = FALSE)
l<-matrix(1,nrow(G),1)

GX<-G %*% x_sim
G2X<-(G %*% G) %*% x_sim
#an identity matrix of appropriate size
I<-(diag(nrow(G))) 
# Inv corresponds to (I - Beta*G))^(-1)   in the reduced form(check equation (5))
# Solve function gives the inverse of a matrix
Inv<-solve(I - beta * G)
```

## Case without fixed effects

```r
%%R

#the instrument vector of size m x 4
S<-matrix(c(l,x_sim,GX,G2X),nrow(x_sim),)

#P is the weighting matrix of size m x n
P<-S %*% solve(t(S) %*% S) %*% t(S)


eps<-matrix(rnorm(n = nrow(G),mean = 0,sd = e_var),nrow(G),1)
y<-alpha * Inv %*% l  + Inv %*% (gamma * I + delta * G) %*% x_sim  + Inv %*% eps
Gy<-G %*% y

#X tilde, size n x 4
X_t<-matrix(c(l,Gy,x_sim,GX),nrow(x_sim),)

#theta 2sls and extracting its parameters
th_2sls <-solve(t(X_t) %*% P %*% X_t) %*% t(X_t) %*% P %*% y
alpha_2sls<-th_2sls[1]
beta_2sls<-th_2sls[2]
gamma_2sls<-th_2sls[3]
delta_2sls<-th_2sls[4]

#Recalculate I with Beta_2sls
I_2sls<-solve((diag(nrow(G)) - beta_2sls * G ))  

#Gy estimated in theta 2sls
gy_2sls<- G %*% I_2sls %*% (alpha_2sls * l   + gamma_2sls * x_sim + GX * delta_2sls)

Z_th<-matrix(c(rep(1,nrow(G)),gy_2sls,x_sim,GX),nrow(x_sim),)

#THETAS
#theta lee
th_lee_1<-solve(t(Z_th) %*% X_t) %*% t(Z_th) %*% y
```

## Case with fixed effects


```r
%%R
IG <-(I - G)
#the instrument vector of size m x 3
S<-matrix(c(IG%*%x_sim,IG%*%GX,IG%*%G2X),nrow(G),)

#P is the weighting matrix of size m x m
P<-S %*% solve(t(S) %*% S) %*% t(S)

y2<- solve(IG) %*% Inv %*% (gamma*I + delta * G) %*% IG %*% x_sim + solve(IG) %*% Inv %*% IG %*% eps

#X tilde, size m x 3
X_t<-matrix(c(G%*%IG%*%y2 ,IG%*%x_sim,IG%*%GX),nrow(x_sim),)


#theta 2sls and extracting its parameters
th_2sls <-solve(t(X_t) %*% P %*% X_t) %*% t(X_t) %*% P %*% IG%*%y2
beta_2sls<-th_2sls[1]
gamma_2sls<-th_2sls[2]
delta_2sls<-th_2sls[3]

#Recalculate I with Beta_2sls
Inv_2sls<-solve(I - beta_2sls * G )
#IGy estimated in theta 2sls
IGy_2sls<-Inv_2sls %*% (gamma_2sls * I + delta_2sls * G) %*% IG %*% x_sim + Inv_2sls %*% IG %*% eps 

IG_Gy_2sls<- G %*% Inv_2sls %*% (IG %*% (x_sim * gamma_2sls + GX* delta_2sls))

Z_th<-matrix(c(IG_Gy_2sls,IG%*%x_sim,IG%*%GX),nrow(G),)


#THETAS
#theta lee
th_lee_2<-solve(t(Z_th) %*% X_t) %*% t(Z_th) %*% IG%*%y2
```

## Storing data in R as CSV, and reading it in Stata


```r
%%R
#storing data and adjacency matrix as CSV, to be readed in Stata later
write.csv(data.frame(x_sim,y,y2),'data.csv',row.names = FALSE)
write.csv(G,'adjacency.csv',row.names = FALSE)
```


```python
%%stata
*reading adjacency matrix, and passing it to Mata
import delimited "adjacency.csv", clear
putmata  G_r=(v*), replace

*reading data
import delimited "data.csv", varnames(1) clear
```
   
    . *reading adjacency matrix, and passing it to Mata
    . import delimited "adjacency.csv", clear
    (encoding automatically selected: ISO-8859-2)
    (100 vars, 100 obs)
    
    . putmata  G_r=(v*), replace
    (1 matrix posted)
    
    . 
    . *reading data
    . import delimited "data.csv", varnames(1) clear
    (encoding automatically selected: ISO-8859-1)
    (3 vars, 100 obs)
    
    . 
    

## Comparing results without fixed effects


```python
%%stata
peer_iv y x_sim, adj(G_r)
```

    -----------------------------------------------------------------------------
    > -
               y | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
    -------------+----------------------------------------------------------------
           _cons |   .7693815   .0861937     8.93   0.000     .5982885    .9404746
             y_p |   .4668116   .0019521   239.13   0.000     .4629367    .4706865
           x_sim |   .0832526   .0174479     4.77   0.000     .0486188    .1178864
         x_sim_p |   .1501907   .0057371    26.18   0.000     .1388026    .1615789
    ------------------------------------------------------------------------------
    


```r
%%R
th_lee_1
```

              [,1]
    [1,] 0.7693815
    [2,] 0.4668116
    [3,] 0.0832526
    [4,] 0.1501907
    

## Comparing results with fixed effects


```python
%%stata
peer_iv y2 x_sim, fixed adj(G_r)
```

    ------------------------------------------------------------------------------
              y2 | Coefficient  Std. err.      t    P>|t|     [95% conf. interval]
    -------------+----------------------------------------------------------------
            y2_p |   .4663327   .0025075   185.97   0.000      .461356    .4713095
           x_sim |   .0841561   .0081916    10.27   0.000      .067898    .1004143
         x_sim_p |   .1500943   .0018714    80.20   0.000       .14638    .1538086
    ------------------------------------------------------------------------------
    


```r
%%R
th_lee_2
```

               [,1]
    [1,] 0.46633273
    [2,] 0.08415615
    [3,] 0.15009431
    

As we can see, in both cases, the `peer_iv` command produces the same estimates as the original package. Furthermore, both packages produce estimates that are very similar to the ones used to generate our data (for instance, we have that $\beta=0.4666$, and in both cases we get coefficients very close to that).

