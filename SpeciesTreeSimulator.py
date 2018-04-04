import AuxiliarFunctions as af
import numpy
import ete3
import random

class SpeciesTreeGenerator():

    def __init__(self, parameters):

        self.parameters = parameters

    def start(self):

        self.whole_tree = ete3.Tree()

        self.lineages_counter = 0
        self.events = list()

        self.active_lineages = set()
        self.inactive_lineages = set()
        self.distances = dict()

        self.active_lineages.add("Root")
        self.distances["Root"] = 0.0
        c1, c2 = self._get_speciated("Root", 0.0)

        return c1, c2

    def run(self):

        self.start()

        speciation = af.obtain_value(self.parameters["SPECIATION"])
        extinction = af.obtain_value(self.parameters["EXTINCTION"])
        stopping_rule = self.parameters["STOPPING_RULE"]
        total_time = self.parameters["TOTAL_TIME"]
        total_lineages = self.parameters["TOTAL_LINEAGES"]
        max_lineages = self.parameters["MAX_LINEAGES"]

        time = 0

        n_lineages_alive = 2

        while True:

            if n_lineages_alive == 0:

                print("All dead")
                success = False
                return success

            print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))
            time_to_next_event = self.get_time_to_next_event(n_lineages_alive, (speciation, extinction))

            if stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                for lineage in self.active_lineages:
                    self.events.append((total_time, "F", lineage))
                success = True
                return success

            elif stopping_rule == 1 and n_lineages_alive == total_lineages:

                self.increase_distances(time_to_next_event)
                for lineage in self.active_lineages:
                    self.events.append((time + time_to_next_event, "F", lineage))
                success = True
                return success

            elif n_lineages_alive >= max_lineages:

                print("Aborting. Max n of lineages attained")
                success = True
                return success

            else:
                # In this case we do the normal the computation

                time += time_to_next_event

                self.increase_distances(time_to_next_event)
                event = self.choose_event(speciation, extinction)
                lineage = random.sample(self.active_lineages, 1)[0]

                if event == "S":
                    self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                elif event == "E":
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1

    def run_b(self):

        # Speciation and extinction rates are branch-wise
        # Each time I create a new lineage, I have to generate new number for its rates
        # We create a dictionary to store the rates

        speciation = af.obtain_value(self.parameters["SPECIATION"])
        extinction = af.obtain_value(self.parameters["EXTINCTION"])

        self.branchwise_rates = dict()
        self.branchwise_rates["Root"] = (speciation, extinction)

        c1, c2 = self.start()

        self.branchwise_rates[c1] = (
            af.obtain_value(self.parameters["SPECIATION"]),
            af.obtain_value(self.parameters["EXTINCTION"]))
        self.branchwise_rates[c2] = (
            af.obtain_value(self.parameters["SPECIATION"]),
            af.obtain_value(self.parameters["EXTINCTION"]))

        stopping_rule = self.parameters["STOPPING_RULE"]
        total_time = self.parameters["TOTAL_TIME"]
        total_lineages = self.parameters["TOTAL_LINEAGES"]
        max_lineages = self.parameters["MAX_LINEAGES"]

        time = 0

        n_lineages_alive = 2


        while True:

            if n_lineages_alive == 0:

                print("All dead")
                success = False
                return success

            print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            time_to_next_event = self.get_time_to_next_event_advanced_modes()

            if stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                for lineage in self.active_lineages:
                    self.events.append((total_time, "F", lineage))
                success = True
                return success

            elif stopping_rule == 1 and n_lineages_alive == total_lineages:

                self.increase_distances(time_to_next_event)
                for lineage in self.active_lineages:
                    self.events.append((time + time_to_next_event, "F", lineage))
                success = True
                return success

            elif n_lineages_alive >= max_lineages:

                print("Aborting. Max n of lineages attained")
                success = True
                return success

            else:
                # In this case we do the normal the computation

                time += time_to_next_event
                self.increase_distances(time_to_next_event)

                # Now we have to choose the lineage doing the event. This will be proportional to the value of the rates
                ###

                active_lineages = list(self.active_lineages)

                lineage = numpy.random.choice(active_lineages, 1, p=af.normalize(
                    [sum((self.branchwise_rates[x][0], self.branchwise_rates[x][1])) for x in active_lineages]))[0]

                myspeciation = self.branchwise_rates[lineage][0]
                myextinction = self.branchwise_rates[lineage][1]

                event = self.choose_event(myspeciation, myextinction)

                if event == "S":
                    c1,c2 = self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                    self.branchwise_rates[c1] = (
                        af.obtain_value(self.parameters["SPECIATION"]),
                        af.obtain_value(self.parameters["EXTINCTION"]))
                    self.branchwise_rates[c2] = (
                        af.obtain_value(self.parameters["SPECIATION"]),
                        af.obtain_value(self.parameters["EXTINCTION"]))

                elif event == "E":
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1


    def run_a(self):

        # Speciation and extinction rates are time autocorrelated
        # Each time I create a new lineage, I have to generate new number for its rates
        # We create a dictionary to store the rates

        speciation = af.obtain_value(self.parameters["SPECIATION"])
        extinction = af.obtain_value(self.parameters["EXTINCTION"])

        self.branchwise_rates = dict()
        self.branchwise_rates["Root"] = (speciation, extinction)

        c1, c2 = self.start()

        self.branchwise_rates[c1] = (
            af.obtain_value(self.parameters["SPECIATION"]),
            af.obtain_value(self.parameters["EXTINCTION"]))
        self.branchwise_rates[c2] = (
            af.obtain_value(self.parameters["SPECIATION"]),
            af.obtain_value(self.parameters["EXTINCTION"]))

        stopping_rule = self.parameters["STOPPING_RULE"]
        total_time = self.parameters["TOTAL_TIME"]
        total_lineages = self.parameters["TOTAL_LINEAGES"]
        max_lineages = self.parameters["MAX_LINEAGES"]

        time = 0

        n_lineages_alive = 2

        while True:

            if n_lineages_alive == 0:

                print("All dead")
                success = False
                return success

            print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            time_to_next_event = self.get_time_to_next_event_advanced_modes()

            if stopping_rule == 0 and time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                for lineage in self.active_lineages:
                    self.events.append((total_time, "F", lineage))
                success = True
                return success

            elif stopping_rule == 1 and n_lineages_alive == total_lineages:

                self.increase_distances(time_to_next_event)
                for lineage in self.active_lineages:
                    self.events.append((time + time_to_next_event, "F", lineage))
                success = True
                return success

            elif n_lineages_alive >= max_lineages:

                print("Aborting. Max n of lineages attained")
                success = True
                return success

            else:
                # In this case we do the normal the computation

                time += time_to_next_event
                self.increase_distances(time_to_next_event)

                # Now we have to choose the lineage doing the event. This will be proportional to the value of the rates
                ###

                active_lineages = list(self.active_lineages)

                lineage = numpy.random.choice(active_lineages, 1, p=af.normalize(
                    [sum((self.branchwise_rates[x][0], self.branchwise_rates[x][1])) for x in active_lineages]))[0]

                myspeciation = self.branchwise_rates[lineage][0]
                myextinction = self.branchwise_rates[lineage][1]

                event = self.choose_event(myspeciation, myextinction)

                if event == "S":
                    c1, c2 = self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                    variance = str(self.distances[lineage] * 0.01)

                    mean_s = str(self.branchwise_rates[lineage][0])
                    mean_e = str(self.branchwise_rates[lineage][1])
                    probability_type_s = self.parameters["SPECIATION"].split(":")[0]
                    probability_type_e = self.parameters["EXTINCTION"].split(":")[0]

                    self.branchwise_rates[c1] = (af.obtain_value(probability_type_s + ":" + mean_s + ";" + variance),
                                                 af.obtain_value(probability_type_e + ":" + mean_e + ";" + variance))

                    print(af.obtain_value(probability_type_s + ":" + mean_s + ";" + variance))

                    mean_s = str(self.branchwise_rates[lineage][0])
                    mean_e = str(self.branchwise_rates[lineage][1])

                    probability_type_s = self.parameters["SPECIATION"].split(":")[0]
                    probability_type_e = self.parameters["EXTINCTION"].split(":")[0]

                    self.branchwise_rates[c2] = (af.obtain_value(probability_type_s + ":" + mean_s + ";" + variance),
                                                 af.obtain_value(probability_type_e + ":" + mean_e + ";" + variance))

                elif event == "E":
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1

    def run_p(self):

        self.start()
        print("Computing tree with fine control of the population size")

        speciation = af.obtain_value(self.parameters["SPECIATION"])
        extinction = af.obtain_value(self.parameters["EXTINCTION"])
        turnover = self.parameters["TURNOVER"]

        time_slices = self.parameters["POPULATION_SIZES"]
        total_time = time_slices[-1][0]
        current_time_slice = 0

        time = 0
        action = 0

        n_lineages_alive = 1

        while True:

            print("Time: %s ; Number of lineages alive: %s" % (str(time), str(n_lineages_alive)))

            goal_N = time_slices[current_time_slice][1]  # The population size we have to attained

            if n_lineages_alive == goal_N:
                time_to_next_event = self.get_time_to_next_event(n_lineages_alive, [turnover])
                action = 0
            elif n_lineages_alive < goal_N:
                time_to_next_event = self.get_time_to_next_event(n_lineages_alive, [speciation])
                action = 1
            elif n_lineages_alive > goal_N:
                time_to_next_event = self.get_time_to_next_event(n_lineages_alive, [extinction])
                action = 2

            if n_lineages_alive == 0:

                self.increase_distances(time_to_next_event)
                print("All dead")
                success = False
                return success

            elif time + time_to_next_event >= total_time:

                self.increase_distances(total_time - time)
                for lineage in self.active_lineages:
                    self.events.append((total_time, "F", lineage))
                success = True
                return success

            else:

                # In this case we do the normal the computation
                time += time_to_next_event
                self.increase_distances(time_to_next_event)

                if time >= time_slices[current_time_slice][0]:
                    current_time_slice +=1

                if action == 0:

                    lineage1, lineage2 = random.sample(self.active_lineages, 2)
                    self._get_speciated(lineage1, time)
                    self._get_extinct(lineage2, time)

                elif action == 1:

                    lineage = random.sample(self.active_lineages, 1)[0]
                    self._get_speciated(lineage, time)
                    n_lineages_alive += 1

                elif action == 2:

                    lineage = random.sample(self.active_lineages, 1)[0]
                    self._get_extinct(lineage, time)
                    n_lineages_alive -= 1

    def increase_distances(self, time):

        for lineage in self.active_lineages:
            self.distances[lineage] += time

    def get_time_to_next_event(self, n, events):

        total = 0.0
        for i in range(n):
            for event in events:
                total += event
        time = numpy.random.exponential(1/total)
        return time

    def get_time_to_next_event_advanced_modes(self):
        # To obtain the time to next event in case that we have different rates per branch
        total = 0.0
        for lineage in self.active_lineages:
            total += self.branchwise_rates[lineage][0]
            total += self.branchwise_rates[lineage][1]
        time = numpy.random.exponential(1 / total)
        return time


    def _get_speciated(self, lineage, time):

        self.lineages_counter += 1
        c1name = "n" + str(self.lineages_counter)
        self.active_lineages.add(c1name)

        self.lineages_counter += 1
        c2name = "n" + str(self.lineages_counter)
        self.active_lineages.add(c2name)

        self.distances[c1name] = 0.0
        self.distances[c2name] = 0.0

        self.inactive_lineages.add(";".join((lineage, c1name, c2name)))
        self.active_lineages.discard(lineage)

        self.events.append((time,"S", ";".join((lineage, c1name, c2name))))  # Store the event

        return c1name, c2name

    def _get_extinct(self, lineage, time):

        self.active_lineages.discard(lineage)
        self.inactive_lineages.add(lineage)
        self.events.append((time, "E", lineage))  # Store the event

    def choose_event(self, speciation, extinction):

        if numpy.random.uniform(0, 1) <= (speciation / (speciation + extinction)):
            return "S"
        else:
            return "E"

    def generate_newick_trees(self):

        def find_descendant(surviving_nodes, node):

            found = 0
            mynode = surviving_nodes[node]["descendant"]

            while found == 0:

                if surviving_nodes[mynode]["state"] == 1:
                    found = 1
                else:
                    mynode = surviving_nodes[mynode]["descendant"]

            return mynode

        # Eric's algorithm

        # First we will iterate the events from the end

        events = self.events

        surviving_nodes = dict()
        times = dict()

        for current_time, event, nodes in events[::-1]:

            if event == "F":

                times[nodes] = float(current_time)
                surviving_nodes[nodes] = {"state": 1, "descendant": "None"}

            elif event == "E":

                times[nodes] = float(current_time)
                surviving_nodes[nodes] = {"state": 0, "descendant": "None"}

            elif event == "S":

                p, c1, c2 = nodes.split(";")

                times[p] = float(current_time)

                if surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == 1:

                    surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + c2}

                elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == 0:

                    surviving_nodes[p] = {"state": 0, "descendant": "None"}

                elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == -1:

                    mynode1 = find_descendant(surviving_nodes, c1)
                    mynode2 = find_descendant(surviving_nodes, c2)

                    surviving_nodes[p] = {"state": 1, "descendant": mynode1 + ";" + mynode2}


                elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == 0:

                    surviving_nodes[p] = {"state": -1, "descendant": c1}

                elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == 1:

                    surviving_nodes[p] = {"state": -1, "descendant": c2}


                elif surviving_nodes[c1]["state"] == 1 and surviving_nodes[c2]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2)
                    surviving_nodes[p] = {"state": 1, "descendant": c1 + ";" + mynode}

                elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == 1:

                    mynode = find_descendant(surviving_nodes, c1)
                    surviving_nodes[p] = {"state": 1, "descendant": mynode + ";" + c2}


                elif surviving_nodes[c1]["state"] == -1 and surviving_nodes[c2]["state"] == 0:

                    mynode = find_descendant(surviving_nodes, c1)
                    surviving_nodes[p] = {"state": -1, "descendant": mynode}

                elif surviving_nodes[c1]["state"] == 0 and surviving_nodes[c2]["state"] == -1:

                    mynode = find_descendant(surviving_nodes, c2)
                    surviving_nodes[p] = {"state": -1, "descendant": mynode}

        extanttree = ete3.Tree()
        wholetree = ete3.Tree()
        eroot = extanttree.get_tree_root()
        eroot.name = ""
        wroot = wholetree.get_tree_root()
        wroot.name = "Root"

        t = (len(events))

        wquick_nodes = dict()
        equick_nodes = dict()

        wquick_nodes["Root"] = wroot

        for i, values in enumerate(events):

            current_time, event, nodes = values

            if event == "S":

                p, c1, c2 = nodes.split(";")

                mynode = wquick_nodes[p]
                myc1 = mynode.add_child()
                myc2 = mynode.add_child()
                myc1.name = c1
                myc2.name = c2
                myc1.dist = times[c1] - times[p]
                myc2.dist = times[c2] - times[p]

                wquick_nodes[c1] = myc1
                wquick_nodes[c2] = myc2

                state = surviving_nodes[p]["state"]

                if state == 1:  # Now the extant tree

                    c1name, c2name = surviving_nodes[p]["descendant"].split(";")

                    if eroot.name == "":
                        eroot.name = p
                        equick_nodes[p] = eroot

                    mynode = equick_nodes[p]

                    myc1 = mynode.add_child()
                    myc2 = mynode.add_child()

                    myc1.name = c1name
                    myc2.name = c2name

                    myc1.dist = times[c1name] - times[p]
                    myc2.dist = times[c2name] - times[p]

                    equick_nodes[c1name] = myc1
                    equick_nodes[c2name] = myc2

        return wholetree.write(format=1), extanttree.write(format=1)

    def write_events_file(self, events_file):

        header = ["TIME","EVENT","NODES"]
        header = "\t".join(map(str, header)) + "\n"

        with open(events_file, "w") as f:
            f.write(header)
            for item in self.events:
                line = "\t".join(map(str,item)) + "\n"
                f.write(line)

    def write_rates(self, rates_file):

        with open(rates_file, "w") as f:

            line = "\t".join(["lineage", "speciation", "extinction"]) + "\n"
            f.write(line)

            for lineage, values in self.branchwise_rates.items():

                speciation, extinction = values

                line = "\t".join(map(str,[lineage, speciation, extinction])) + "\n"
                f.write(line)

