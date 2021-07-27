from polymer import Polymer
import math
import csv
import random

def rr_extracter(rr_df, monomers):
    listy = [[rr_df[a][b] if a != b else 'same' for a in monomers] for b in monomers]
    for i in listy:
        while 'same' in i:
            i.remove('same')
    return listy

def run_simulation(N_MONs, N_CHAINs, MRs, RRs, avgDP, conv, CTP, PRUNE_OLIGOMERS):
    
    
    def getMonomerAmounts(MRs, N_MONs, MP):
        "***Calculates the number of each monomer given poolSize and monomerRatios***"
        monomerAmounts = []
        sumMonomerRatios = sum(MRs)

        for i in range(N_MONs):
            monomerAmounts.append(int(MRs[i]/sumMonomerRatios*MP))

        return monomerAmounts

    
    def testAssertions(numPolymers, avgDP, MRs):
        "***Makes sure inputs are valid***"
        inputsValid = True

        if MP < avgDP/100:
            print("Error", "Monomer Pool Size is too small!")
            inputsValid = False

        elif numPolymers < 1:
            print("Error", "DP or Monomer Pool Size is too small!")
            inputsValid = False

        elif int(avgDP) <= 0:
            print("Error", "DP or Conversion cannot be zero!")
            inputsValid = False

        #elif min(MRs) == 0:
        #    print("Error", "Monomer Ratios cannot be zero!")
        #    inputsValid = False

        return inputsValid
    

    def getRateConstants(N_MONs, RRs):
        "***Get rate constants from reactivity ratio inputs***"
        rateConstantList = []

        for i in range(N_MONs):
            rateConstantList.append([])
            jIndex = 0

            for j in range(N_MONs):

                if i == j:
                    rateConstantList[i].append(1.0)

                else:
                    rrValue = RRs[i][jIndex]
                    if rrValue == 0:
                        rateConstantList[i].append(math.inf)
                    else:
                        rateConstantList[i].append(1/rrValue)
                    jIndex += 1

        return rateConstantList


    def weighted_choice(choices):
        "***Takes in a list of tuple lists with item and weight, and returns a random item based on weight***"
        total = sum(w for c, w in choices)
        r = random.uniform(0, total)
        upto = 0
        for c, w in choices:
            if upto + w >= r:
                return c
            upto += w
        assert False, "Shouldn't get here"

        
    def exportPolymerArray(polymerArray, N_MONs, N_CHAINs, MRs, RRs, avgDP, conv, CTP, PRUNE_OLIGOMERS):
        csvData = []
        for polymer in polymerArray:
            polymer = polymer.asArray()
            if len(polymer) > PRUNE_OLIGOMERS:
                csvData.append(polymer)
         
        file = 'outputs/' + 'NM' + str(N_MONs) + '_MR' + str(MRs) + '_NC' + str(N_CHAINs) + '_DP' + str(avgDP) + '_conv' + str(int(conv*100)) + '_CTP' + str(int(CTP*100)) + '_FILT' + str(PRUNE_OLIGOMERS) + ".csv"
        if file:
            with open(file, 'w') as csvFile:
                writer = csv.writer(csvFile)
                writer.writerows(csvData)
            csvFile.close()
            print("Simulation done! Sequence file saved @", file)
    
    
    "***Calculate monomer pool size***"
    MP = int( N_CHAINs * (avgDP/conv) )
    
    #initiate a Polymer Array, which will hold all of the Polymer Objects to be simulated and created
    polymerArray = []

    #each monomer amount is calculated based off its ratio and the pool size
    monomerAmounts = getMonomerAmounts(MRs, N_MONs, MP)

    #Pool size is adjusted to account for rounding errors while calculating monmer amounts
    adjustedPoolSize = sum(monomerAmounts)

    numPolymers = int(adjustedPoolSize / (avgDP/conv))


    "***Check if user inputs are valid***"
    inputsValid = testAssertions(numPolymers, avgDP, MRs)
    if not inputsValid:
        return

    
    # get rate constants from RRs
    rateConstantList = getRateConstants(N_MONs, RRs)

    "---------------------------------------------------------------------------------------------------------------------------------"
    "***INITIATION STEP***"
    "---------------------------------------------------------------------------------------------------------------------------------"

    "In this step, all polymers will be initiated with one starting monomer. The algorithm continually intiates"
    "a chain until the chain count reaches numPolymers"

    for i in range(numPolymers):  

        "In the case where the user fixes all chains to start with a certain monomer, we simply initiate all chains"
        "with the monomer of choice."

        #case for setting first monomer
        "Otherwise, we will choose the initial monomer randomly based on the the mole fractions 'f' of each monomer"
        "in the feed solution. For 2- and 3- monomer systems, the instantaneous form of the Mayo Lewis Equation is "
        "used to determine the probabilities of each monomer initiating the chain."

        "Initiate a variable 'choices' that keeps track of the weight probabilty assigned to each monomer."
        #Example:

        #>>>print(choices)
        #[[1, 1.5], [2, 0.5]]

        #This means that monomer 1 has a weight of 1.5 assigned to it and monomer 2 has a weight of 0.5 assigned to it.
        #Thus, monomer 1 will initiate 3x more often than monomer 2 will.

        choices = []


        if N_MONs == 2:
            "Case for 2-monomer system: use the Mayo Lewis Equation"
            f1 = monomerAmounts[0]
            f2 = monomerAmounts[1]
            r1 = RRs[0][1]
            r2 = RRs[1][0]
            weight = (r1*f1**2 + f1*f2) / (r1*f1**2 + 2*f1*f2 + r2*f2**2)
            choices.append([1, weight])
            choices.append([2, 1 - weight])

        elif N_MONs == 3:
            "Case for 3-monomer system: use an altered Mayo Lewis Equation"
            m1 = monomerAmounts[0]
            m2 = monomerAmounts[1]
            m3 = monomerAmounts[2]
            F = m1 + m2 + m3
            f1 = m1/F
            f2 = m2/F
            f3 = m3/F
            r11 = rateConstantList[0][0]
            r12 = rateConstantList[0][1]
            r13 = rateConstantList[0][2]
            r21 = rateConstantList[1][0]
            r22 = rateConstantList[1][1]
            r23 = rateConstantList[1][2]
            r31 = rateConstantList[2][0]
            r32 = rateConstantList[2][1]
            r33 = rateConstantList[2][2]
            R1 = r11 + r12 + r13
            R2 = r21 + r22 + r23
            R3 = r31 + r32 + r33
            a = f1*r11*f1/(r11*f1+r12*f2+r13*f3) + f2*r21*f1/(r21*f1+r22*f2+r23*f3) + f3*r31*f1/(r31*f1+r32*f2+r33*f3)
            b = f1*r12*f2/(r11*f1+r12*f2+r13*f3) + f2*r22*f2/(r21*f1+r22*f2+r23*f3) + f3*r32*f2/(r31*f1+r32*f2+r33*f3)
            c = 1 - a - b
            choices.append([1,a])
            choices.append([2,b])
            choices.append([3,c])

        else:
            "Case for any monomer system > 3: initial solely based on feed mole ratios 'f'"
            "Cycle through each monomer, finding its initial amount and using that value as the weight."

            for i in range(N_MONs):
                #weight chance of monomer initation: (amount of starting monomer)
                weight = monomerAmounts[i]
                #Adds a two element list to choices containing monomer and weight
                choices.append([i+1, weight])

                
        "Randomly choose a monomer to initiate the polymer chain by using a weighted random selector."
        "Monomer with higher relative weights will be chosen more often. The 'weighted_choices' function"
        "takes in the 'choices' variable which has relevant weights for each monomer and runs a weighted"
        "random selection on it."
        try:
            startingMonomer = weighted_choice(choices)

        #If all weights are zero due to coefficients all being zero, then sort by relative amounts of monomer instead
        except AssertionError:
            monomerID = 1
            choices = []
            while monomerID <= N_MONs:
                choices.append([monomerID, monomerAmounts[monomerID - 1]])
                monomerID += 1
            startingMonomer = weighted_choice(choices)

        #Starts a new polymer with startingMonomer, represented by an array, 
        #and adds that array to polymerArray
        "A polymer is represented as a Python Class, as seen in Polymer.py. Here, we inititate an instance of the Polymer object, "
        "and append that Polymer to a list containing all of the polymers"
        polymer = Polymer([startingMonomer], N_MONs, rateConstantList)

        polymerArray.append(polymer)

        "Remove the monomer that was used to intitate the polymer in order to accurately update the monomer pool counts."
        monomerAmounts[startingMonomer - 1] -= 1


    "------------------------------------------------------------------------------------------------------------------------"
    "***PROPAGATION STEP***"
    "------------------------------------------------------------------------------------------------------------------------"

    #Calculating how many monomers the reaction will use, taking into account initiated monomers
    monomersUsed = 1*numPolymers

    monomers_to_consume = int(conv * adjustedPoolSize) - monomersUsed

    "***Propogation for standard Mayo-Lewis Case***"

    "The 'for' loop will add monomers to the chain until we reach we use up a percentage of the total starting monomers"
    "defined by user input into the 'Percent Conversion' box. Setting 'Percent Conversion' to 100% will continue to add "
    "monomers to growing chains until there are no remaining monomers, simulating a living polymerization."

    i = 0;
    while i < monomers_to_consume:

        "Randomly choose a polymer chain to grow."
        polymer = random.choice(polymerArray)

        "For the polymer selected, iterate through all possible monomers which can be added. For each monomer, calculate"
        "the weight chance of the monomer to be added to the chain , defined as the product of the relevant rate "
        "constant 'k' times the number of monomers left unreacted 'f'"

        #A varaible to keep track of number of monomers to add. This is usually 1; however if the Chain Transfer % is less than 100, there
        #is a chance that more than one monomer is added.
        num_monomers_to_add = 1

        #percentage that an extra monomer will be added
        fudge_factor = 1 - CTP

        #calculate how many monomers to add
        while True:
            #Randomly determine if another monomer will be added
            if random.random() >= fudge_factor:
                break
            #If number of monomers added exceeds total monomers left to consume, stop
            elif num_monomers_to_add >= monomers_to_consume - num_monomers_to_add - i:
                break
            else:
                num_monomers_to_add += 1

        for j in range(num_monomers_to_add):
            #A variable keeping track of monomer choices and weights
            choices = []

            for monomerID in range(1, N_MONs + 1):
                #Retrieveing coefficient based on previous and current monomer

                #retrieve the relevant rate constant k
                k = polymer.rateConstant(monomerID)

                # weight chance calulations for monomer attaching: coefficient*(amount of monomer remaining)
                chance = monomerAmounts[monomerID - 1] * k
                #Adds a two element list to choices containing monomer and weight
                choices.append([monomerID, chance])

            "Using weighted_choice, select next monomer to be appended to the growing chain"
            try:
                nextMonomer = weighted_choice(choices)

            #If all weights are zero due to coefficients all being zero, then sort by relative amounts of monomer instead
            except AssertionError:
                monomerID = 1
                choices = []
                while monomerID <= N_MONs:
                    choices.append([monomerID, monomerAmounts[monomerID - 1]])
                    monomerID += 1
                nextMonomer = weighted_choice(choices)

            "Attach next monomer to polymer chain"
            polymer.append(nextMonomer)

            "Remove the selected monomer from the pool to reflect the monomer being used up in the reaction"
            "If Hold Composition is checked, then the simulation will not use up any monomers, and end when the expected number of"
            "mononomers are consumed instead"
            monomerAmounts[nextMonomer - 1] -= 1

            "Increment counter for monomers used"
            i += 1

    exportPolymerArray(polymerArray, N_MONs, N_CHAINs, MRs, RRs, avgDP, conv, CTP, PRUNE_OLIGOMERS)

    "***Update Global Variables***"
    #lambdaValue = analysis.calculate_theta(self) - do we care about the lambdavalue?

