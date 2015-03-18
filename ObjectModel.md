# Introduction #

In order for the code to be more reusable and testable, I will write it in an object-oriented manner. This page lists the Objects I will implement, their methods and proerties, and the way they interact among each other.

---

# Objects #

## Network ##
This object contains the conectivity data for a single network. It takes care of the logic evaluation, wiring cost evaluation, modularity calculation, etc. It has functions that make create point mutations, and a random network constructor.
## Population ##
This is a collection of networks. It can apply the functions of the single networks in bulk and store the results. It has a subscript operator in order to address individual networks.
## Experiment ##
This is the highest level object. it contains the genetic algorithm, instanciates the populations and handles data output.