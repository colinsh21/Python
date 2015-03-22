# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 15:19:42 2015

@author: colinsh
"""

# Imports
import copy
import networkx as nx
import numpy as np
import pandas
import seaborn; seaborn.set()

class Model(object):
    def __init__(self,graph,initial_infected=[[366],[1000]],reproduction_number=.65,effective_contacts=15,prob_of_flying_day=.0005,passengers_per_flight=150):
        #set pararms
        self.g=graph
        self.initial_infected=initial_infected
        self.reproduction_number=reproduction_number
        self.prob_of_flying_day=prob_of_flying_day
        self.prob_of_flying=1-(1-self.prob_of_flying_day)**30
        self.effective_contacts=effective_contacts
        self.transmission_rate=self.reproduction_number/self.effective_contacts
        self.passengers_per_flight=passengers_per_flight
        
        #set state variables
        self.t=0
        self.history_infected=[]
        self.history_cities_infected=[]
        self.history_infected_interactions=[]
        self.history_space=[]
        self.history_infected_interactions=[]
        self.num_interactions=0
        self.history_num_interactions=[]
        self.time_series_data=[]
        #self call setup
        self.setup_space()
        
    def setup_space(self):
        self.num_infected=0
        self.num_hosts_infected=0
        for node_id in self.g.nodes():
            if node_id in self.initial_infected[0]:
                self.g.node[node_id]['infected']=self.initial_infected[1][0]
                self.g.node[node_id]["state"] = "I"
                self.num_hosts_infected+=1
                self.num_infected+=self.g.node[node_id]['infected']
            else:
                self.g.node[node_id]['infected']=0
                self.g.node[node_id]["state"] = "S"
        
        # Track the latest step
        self.history_space.append(copy.deepcopy(self.g))
        self.history_infected.append(self.num_infected)
        self.history_cities_infected.append(self.num_hosts_infected)
        self.num_interactions=0
        self.time_series_data.append((self.num_infected,self.num_hosts_infected,self.num_interactions))
        
        
    def step(self):
        #self.g_step=self.history_space[-1]
        self.g_step=self.g
        self.num_interactions=0
        # Iterate over I and infect any S neighbors
        for self.node_id in self.g.nodes():
            if self.g.node[self.node_id]["state"] == "I":

                # Decide how many infected are traveling
                self.infected_traveling=np.random.binomial(self.g.node[self.node_id]['infected'],self.prob_of_flying)

                #check
                if self.infected_traveling>self.g.node[self.node_id]['infected']:
                    self.infected_traveling=self.g.node[self.node_id]['infected']
                    
                # Get travel options
                self.neighbors = self.g.neighbors(self.node_id)
                #print('neighbors',neighbors)
                self.passenger_list=[]
                for self.neighbor_id in self.neighbors:
                    self.passenger_list.append(self.g.edge[self.node_id][self.neighbor_id]['weight'])

                #print('passenger_list',passenger_list)
                # Create passenger distribution
                self.percent_list = [float(i)/sum(self.passenger_list) for i in self.passenger_list]

                #passenger_dist=np.cumsum(percent_list)
                #print('passenger_dist',percent_list)

                # Create travel distribution
                self.travel_dist=np.random.multinomial(self.infected_traveling, self.percent_list)
                #print('travel_dist',travel_dist)

                # Convert number of travels to infection
                self.infected_travelers=[]
                for i in xrange(len(self.travel_dist)):
                    self.moving_passengers=self.travel_dist[i]
                    self.max_passengers=self.g.edge[self.node_id][self.neighbors[i]]['weight']
                    self.moving_passengers*=self.transmission_rate*self.passengers_per_flight*np.random.normal(1,.15,1)
                    self.moving_passengers=int(self.moving_passengers)
                    if self.moving_passengers>self.max_passengers:
                        self.moving_passengers=self.max_passengers
                        
                    self.infected_travelers.append(self.moving_passengers)    
                #print self.infected_travelers
                
                # Node update
                # Calculate in host spread
                #print self.node_id,'initial',self.g.node[self.node_id]['infected']
                self.new_infected=int(self.g.node[self.node_id]['infected']*self.reproduction_number)
                #print self.node_id,'spread',self.new_infected
                # Calculate in host death
                self.removed_infected=self.g.node[self.node_id]['infected']
                #print self.node_id,'removed',self.removed_infected
                #print self.node_id,'flying',sum(self.infected_travelers)
                # Spread to neighbors
                for i in xrange(len(self.neighbors)):
                    self.neighbor_id=self.neighbors[i]

                    # Move Infected
                    self.num_infected_at_neighbor=self.g.node[self.neighbor_id]['infected']
                    self.g_step.node[self.neighbor_id]['infected']=self.num_infected_at_neighbor+self.infected_travelers[i]
                    if self.infected_travelers[i]>0:
                        self.history_infected_interactions.append((self.t,self.node_id,self.neighbor_id,self.travel_dist[i],self.infected_travelers[i]))
                        self.num_interactions+=1
                    #self.g_step.node[self.node_id]['infected']=self.g.node[self.node_id]['infected']-self.infected_travelers[i]

                # Sum new infection change
                self.new_infected_total=int(max(self.g.node[self.node_id]['infected']+self.new_infected-self.removed_infected,0))
                self.g_step.node[self.node_id]['infected']=int(self.new_infected_total)
                #print self.node_id,'new infected total',self.new_infected_total
                #print self.node_id,'new infected total graph',self.g_step.node[self.node_id]['infected']
                
        # Update graph
        self.num_infected=0
        self.num_hosts_infected=0
        for self.node_id in self.g_step.nodes():
            if self.g_step.node[self.node_id]['infected']>0:
                #print self.node_id,'END pre-count number infected',self.g_step.node[self.node_id]['infected']
                self.g_step.node[self.node_id]["state"] = "I"
                self.num_hosts_infected+=1
                self.num_infected+=self.g_step.node[self.node_id]['infected']
                #print self.node_id,'number infected',self.g_step.node[self.node_id]['infected']
            else:
                self.g_step.node[self.node_id]["state"] = "S"

        # Resolve g_step to g ***Not Resolving***
        self.g=self.g_step

        # Track the latest step
        self.history_space.append(copy.deepcopy(self.g))
        self.t+=1
        self.history_infected.append(self.num_infected)
        self.history_cities_infected.append(self.num_hosts_infected)
        self.history_num_interactions.append(self.num_interactions)
        self.time_series_data.append((self.num_infected,self.num_hosts_infected,self.num_interactions))
        
