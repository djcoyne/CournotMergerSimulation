# CournotMergerSimulation
A simulation program for Mergers in Cournot Markets

This Python program simulates a merger in Cournot Markets. It allows for a user-defined number of competitors and distinct geographic markets. (Functionality for jointly-defined geographic and product markets will be added in a future update.)

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
         
      
    cfirm - These are firm objects that are passed into a CournotMerge object
  
      Functions:
        __init__(self, name="", q=[0], markets=[""])
          returns: None
          
          Initializes the cfirm object
        
        describe(self)
          returns: string
          
          Describes production for the firm by market
        
        merge(self, partner)
          returns: cfirm object
          
          Merges self with partner by combining names and summing quantities by market
         
         
sample.py - Contains an imaginary merger of avocado retailers that compete in 3 geographic markets. Demonstrates how a user-made  script can use the CournotMerge program
    
