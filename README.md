# unexpected_improv_of_speedy_algo

This unfinished investigation found that utilizing state decompositions in the greedy algorithm proposed in ["Speedy Contraction of ZX Diagrams with Triangles via Stabiliser Decompositions"](https://arxiv.org/abs/2307.01803) is only beneficial
in certain cases, and might hurt the performance in others (where instead dynamic decompositions are used).

In particular, this [log file](https://github.com/Tix3Dev/unexpected_improv_of_speedy_algo/blob/main/zxbarren-private/result.txt) seems to suggest that the use is justified for the application to quantum machine learning ansätze, as on average
it performs much better than the version that does not use state decompositions (in certain cases it is slightly worse, but in most cases it is either as good or significantly better).

On the other hand, when using the version without state decompositions for the tested randomly generated multi-control Toffoli gate dense quantum circuits, we can observe reliable improvements over the original implementation:

<img src="https://github.com/Tix3Dev/unexpected_improv_of_speedy_algo/blob/main/average_performance_ratio_result.png">

<img src="https://github.com/Tix3Dev/unexpected_improv_of_speedy_algo/blob/main/percentage_of_improvements_result.png">

Although the exact reasons for this behavior would need to be investigated, one possible explanation could be the following. When using state decompositions, the resulting terms contain newly added states, consisting of X-spiders, and more importantly, Z-spiders.
In contrast, dynamic decompositions only contain X-spiders. They have the pleasant property of being able to turn star edges into Clifford states, which turns out to be a common pattern. Although technically possible, this occurs a lot less for Z-spiders.
A potential reason for why this could positively affect multi-control Toffoli gates is that those circuits tend to be highly-interconnected. Therefore, potential simplications resulting from X-spiders could outweigh a slight advantage of state decompositions
with regard to their scaling. Conversly, this might not be the case for the less interconnected quantum machine learning diagrams.
