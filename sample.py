'''
This program provides an example of how to use the cournotMerge simulation.
It explores an imaginary merger of avocado retailers in 3 markets.
'''

import CournotMerge as cm

# Define the markets (use the same order for elasticities, prices, and firm quantities below)
markets = ["San Diego", "Seattle", "DC"]

# Define the market elasticities (order same as markets list)
elasticity = [0.75, 1.6, 3.1]

# Define the market prices
price = [0.66, 1.20, 1.90]

# Define the Firms as a dictionary (Name, quantities, markets)
firms={}
firms[1] = cm.firm("Merchant Rick's", [45, 25, 15], markets)
firms[2] = cm.firm("Dubertson's", [40, 18, 6], markets)
firms[3] = cm.firm("Half Foods", [55, 45, 19], markets)
firms[4] = cm.firm("Seedlings", [11, 7, 4], markets)
firms[5] = cm.firm("Krugel", [38, 23, 10], markets)

# Identify the keys for the two parties that are merging
merge_keys = ['MR','HF']

# Identify any MC efficiencies for the merging parties (efficiency% reduction in MC)
efficiency = 5

# Call the merger program
CM = cm.cournotMerge(F=firms, markets=markets, e=elasticity, p=price, MF=merge_keys, E=efficiency)

# Describe the merger simulation
print(CM.describe())

# Summarize the merger sim
CM.summarize()
