# Variational-Self-Organizing-Maps
A c++ implementation of the original SOM algorithm by Teuvo Kohonen, but with a variational twist to help with anomaly detection. Also features a SOM version of the variational autoencoder

## Variational modeling
In the original Self-Organizing Maps algorithm, every model vector is estimating the mean of a clustered subset of the data:

$$M_i \approx \frac{1}{N}\sum^S v_n$$

where $M_i$ is a model vector, $S$ is the subset, $v_n$ is a datapoint that belong to that subset and $N$ is the amount of such datapoints.

The VSOM (Variational Self-Organizing Maps) also tries to model $M_i$ as a multivariate gaussian distribution with a diagonal covariance matrix, i.e. all dimensions are independent. So, while also estimating the mean as in the original SOM algorithm, the VSOM also estimates:

$$\sigma_i^2 \approx \frac{1}{N-1}\sum^S (v_n - M_n)^2$$

where $\sigma_i^2$ is the model vector $M_i$'s variance.


## Transformations
The VSOM library uses a ``Transformation`` object that is injected into the algorithm in order to transform the data and specify how operations are to be applied to it. At the moment, two non-standard tranformations are supplied.

### Median estimator
The SOM standard update algorithm calculates a new $M_i^{t+1}$ in the following way:

$$M_i^{t+1} = M_i^t + \alpha*(v_i - M_i^t)$$

where $M_i^t$ is the initial model vector, $v_i$ is an observed vector and $\alpha$ is some (small) scalar step size. That will converge to the mean value of the data. The Median estimator transform, on the other hand, uses the following update equation:

$$M_i^{t+1} = M_i^t + \alpha*sign(v_i - M_i^t)$$

where $sign(x)$ equals $-1$ if $x$ is negative, $1$ if $x$ is positive and $0$ otherwise. This equation will make $M_i$ converge towards the median of the clustered data.

### Combinatorial linear regression
This transform models and clusters linear *relationships* between the dataset's dimensions. Given a dataset with $J$ dimensions, then 

$$v_i^j$$

is the value of dimension $j$ of observation $i$. $v_i^j$ has $J-1$ relations to the other dimensions than $j$. For the whole observed vector $v_i$, there are in total

$$\frac{J(J-1)}{2}$$

number of *unique* relations between different dimensions. Each of these relations is modeled with a linear equation:

$$v_i^{j*} = Av_i^k + B$$

where $v_i^k$ is the value of dimension $k$ of an observation $i$, $v_i^{j*}$ is the *modeled* value of dimension $j$ of the same observation $i$ and $A$ and $B$ are two parameters to be fit to the training data. The total number of parameters is thus:

$$N_p = \frac{J(J-1)}{2}*2 = J(J-1)$$

On each update, the linear regression learning objective is followed:

$$J = \sum_{S,R}(Av_i^{r_1}+B - v_i^{r_2})^2$$

where $S$ is a subset, or cluster that observation $v_i$ is assigned to, $R$ is out *relation* space, $v_i^{r_1}$ is an observed dimension of $v_i$ and $v_i^{r_2}$ is the dependent dimension. This boils down to the following update equation:

$$M_i^{A,t+1} = M_i^t + 2\alpha(Av_i^{r_1} + B - v_y^{r_2})v_i^{r_1}$$
$$M_i^{B,t+1} = M_i^t + 2\alpha(Av_i^{r_1} + B - v_y^{r_2})$$

where $M_i^A$ is the proportionality parameter of the model vector and $M_i^B$ is the bias parameter of the model vector.

Since each model vector represent the internal *relations* of a subset of the training set, clustering with the Combinatorial linear regression transform occurs in the *relation* space. This way, other insights can be achived from the data than with the standard SOM algorithm.

# Installation
## Dependencies
### Eigen
'''sudo apt install libeigen3-dev

cd /usr/include
sudo ln -sf eigen3/Eigen Eigen
sudo ln -sf eigen3/unsupported unsupported'''