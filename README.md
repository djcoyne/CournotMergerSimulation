# CournotMergerSimulation
A simulation program for Mergers in Cournot Markets

This Python program simulates a merger of 2 firms in an imperfectly competitive industry. It generates predicted effects of a merger/acquisition by leveraging the following
assumptions:

a. Cournot competition between user-defined firms
b. Linear demand, potentially in multiple geographic markets
c. quadratic costs

(Note: The program also currently assumes that production is at a single facility and that firms compete over one specific product market. Functionality for jointly-defined geographic and product markets, as well as support for multiple production facilities will be added in a future update.)

Each merger simulation requires the following inputs:
1. Matrix of production quantities/capacities at the firm level
2. Market demand elasticity
3. Pre-merger price
4. Any merger-specific marginal cost efficiency

# Files included:
CournotMerge.py - Contains the classes and functions required to run the merger simulation
  
  Classes:
    
    CournotMerge - These are merger objects that describe a specific merger scenario
    
      Functions:
      
        __init__(self, F={}, markets=[], e=[], p=[], MF = [0,1], E=0)
          returns: None
          
          Initializes the CournotMerge object and runs the simulation program
        
        describe(self)
          returns: string
          
          Describes the merger the object analyzes 
          
         summarize(self)
          returns: None 
          
          Summarizes the merger by displaying a host of useful information and analysis
          
         production(self)
          returns: string
          
          Displays the post-merger production and percent change by market
          
         prices(self)
          returns: string
          
          Displays the post-merger prices and percent change by market
         
         concentration(self)
          returns: string
          
          Displays the post-merger HHI and delta by market
          
         profits(self)
          returns: string
          
          Displays the post-merger profits and percent change by firm
          
         demand(self)
          returns: string
          
          Displays the inverse demand functions for each market
         
         harm(self)
          returns: string
          
          Displays the estimated harm to consumer surplus as a result of the merger by market
         
      
    firm - These are firm objects that are passed into a CournotMerge object
  
      Functions:
        __init__(self, name="", q=[0], markets=[""])
          returns: None
          
          Initializes the cfirm object
        
        describe(self)
          returns: string
          
          Describes production for the firm by market
        
        merge(self, partner)
          returns: firm object
          
          Creates a new firm by merging self with partner (combines names and sums production quantities by market)
         
         
sample.py - Contains an imaginary merger of avocado retailers that compete in 3 geographic markets. Demonstrates how a user-made  script can use the CournotMerge program
    
