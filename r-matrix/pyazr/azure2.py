import os
import time

import numpy as np

from .client import client
from .server import server

from concurrent.futures import ThreadPoolExecutor

class azure2:

    servers = []

    def __init__( self, file, nprocs = 1, port = 20000 ):
        self.PORT = port
        self.file= file
        self.clients = []
        self.spawn( nprocs )
        self.configure( )

    def __del__( self ):
        pass

    def spawn( self, nprocs ):
        # Start all servers in parallel
        def _start_server(i):
            return server( self.PORT + i, self.file )
        with ThreadPoolExecutor(max_workers=nprocs) as ex:
            self.servers = list(ex.map(_start_server, range(nprocs)))
        time.sleep( 1 )

        # Connect all clients in parallel
        def _start_client(i):
            return client( port=self.PORT + i )
        with ThreadPoolExecutor(max_workers=nprocs) as ex:
            self.clients = list(ex.map(_start_client, range(nprocs)))
        time.sleep( 1 )

        # Initialize all clients in parallel
        def _initialize(c):
            c.communicate( "INITIALIZE", [0] )
        with ThreadPoolExecutor(max_workers=nprocs) as ex:
            list(ex.map(_initialize, self.clients))
        time.sleep( 1 )

    def configure( self ):
        self.nsegments = int(self.clients[0].communicate( "UPDATE_DATA", [0] )[0])
        self.energies = [ self.clients[0].communicate( "GET_DATA_ENERGIES", [i] ) for i in range( self.nsegments ) ]
        self.excitation_energies = [ self.clients[0].communicate( "GET_DATA_EXCITATION_ENERGY", [i] ) for i in range( self.nsegments ) ]
        self.angles = [ self.clients[0].communicate( "GET_DATA_ANGLES", [i] ) for i in range( self.nsegments ) ]
        self.cross = [ self.clients[0].communicate( "GET_DATA_SEGMENTS", [i] ) for i in range( self.nsegments ) ]
        self.cross_err = [ self.clients[0].communicate( "GET_DATA_SEGMENTS_ERRORS", [i] ) for i in range( self.nsegments ) ]
        self.conv = [ self.clients[0].communicate( "GET_DATA_CONV", [i] ) for i in range( self.nsegments ) ]
        self.sfactor = [ self.cross[i] * self.conv[i] for i in range( self.nsegments ) ]
        self.sfactor_err = [ self.cross_err[i] * self.conv[i] for i in range( self.nsegments ) ]
        self.params = self.clients[0].communicate( "GET_PARAMS", [0] )
        self.params_rwa = self.clients[0].communicate( "GET_PARAMS_RWA", [0] )
        self.fixed_params = self.clients[0].communicate( "GET_PARAMS_FIXED", [0] )

    def norm_indices( self, proc = 0 ):
        indices = self.clients[proc].communicate( "GET_NORM_INDICES", [0] )
        return indices
    
    def shift_indices( self, proc = 0 ):
        indices = self.clients[proc].communicate( "GET_SHIFT_INDICES", [0] )
        return indices

    def update_rwa_params_from_sav( self, dir = '' ):
        if dir == '': 
            all_rwa_params = np.loadtxt('output/param.sav', usecols=(1,))
            all_rwa_names = np.loadtxt('output/param.sav', usecols=(0,), dtype=str)
        else: 
            all_rwa_params = np.loadtxt(dir, usecols=(1,))
            all_rwa_names = np.loadtxt(dir, usecols=(0,), dtype=str)
        self.params_rwa = []
        self.params_names =[]
        for i in range(len(all_rwa_params)):
            if self.fixed_params[i]:
                continue
            else:
                self.params_rwa.append(all_rwa_params[i])
                self.params_names.append(all_rwa_names[i])
    
    def update_sav_from_rwa_params( self, best ):
        params_full = []
        with open('output/param.sav', 'r') as f:
            lines = f.readlines()
            for i in range(len(lines)):
                l = lines[i].split()
                param = [l[0], float(l[1]), float(l[2])]
                params_full.append(param)
        
        # Now update the params_full with best for non-fixed params
        idx = 0
        for i in range(len(params_full)):
            if self.fixed_params[i]: continue
            else:
                params_full[i][1] = best[idx]
                idx += 1

        # Write back to param.sav
        with open('output/param.sav.new', 'w') as f:
            for param in params_full:
                f.write(f"{param[0]} {param[1]} {param[2]}\n")

    def transform_rwa( self, params, proc = 0 ):
        params_physical = self.clients[proc].communicate( "TRANSFORM_RWA", params )
        return params_physical

    def transform_physical( self, params, proc = 0 ):
        params_physical = self.clients[proc].communicate( "TRANSFORM_PHYSICAL", params )
        return params_physical
    
    def calculate_excitation_energy( self, params, proc = 0 ):
        nsegments = int(self.clients[proc].communicate( "UPDATE_SEGMENTS", params ))
        segments = [ self.clients[proc].communicate( "GET_EXCITATION_ENERGY", [i] ) for i in range( nsegments ) ]
        return segments

    def calculate_angles( self, params, proc = 0 ):
        angles = self.clients[proc].communicate( "GET_DATA_ANGLES", params )
        return angles

    def transform_all_rwa( self, params, proc = 0 ):
        params_physical = self.clients[proc].communicate( "TRANSFORM_ALL_RWA", params )
        return params_physical
    
    def calculate_chi2_rwa( self, params, proc = 0 ):
        chi2 = self.clients[proc].communicate( "CALCULATE_CHI2_RWA", params )
        return chi2.tolist()

    def calculate_chi2( self, params, proc = 0 ):
        chi2 = self.clients[proc].communicate( "CALCULATE_CHI2", params )
        return chi2.tolist()

    def calculate( self, params, proc = 0 ):
        nsegments = int(self.clients[proc].communicate( "UPDATE_SEGMENTS", params ))
        segments = [ self.clients[proc].communicate( "GET_CALCULATED_SEGMENT", [i] ) for i in range( nsegments ) ]
        return segments
    
    def calculate_rwa( self, params, proc = 0 ):
        nsegments = int(self.clients[proc].communicate( "UPDATE_SEGMENTS_RWA", params ))
        segments = [ self.clients[proc].communicate( "GET_CALCULATED_SEGMENT", [i] ) for i in range( nsegments ) ]
        return segments
    
    def calculate_all_rwa( self, params, proc = 0 ):
        nsegments = int(self.clients[proc].communicate( "UPDATE_SEGMENTS_ALL_RWA", params ))
        segments = [ self.clients[proc].communicate( "GET_CALCULATED_SEGMENT", [i] ) for i in range( nsegments ) ]
        return segments
    
    def calculate_energies( self, params, proc = 0 ):
        nsegments = int(self.clients[proc].communicate( "UPDATE_SEGMENTS", params ))
        segments = [ self.clients[proc].communicate( "GET_CALCULATED_ENERGIES", [i] ) for i in range( nsegments ) ]
        return segments
    
    def calculate_sfactor( self, params, proc = 0 ):
        nsegments = int(self.clients[proc].communicate( "UPDATE_SEGMENTS", params ))
        segments = [ self.clients[proc].communicate( "GET_CALCULATED_SEGMENT", [i] ) for i in range( nsegments ) ]
        conv = [ self.clients[proc].communicate( "GET_CALCULATED_CONV", [i] ) for i in range( nsegments ) ]
        for i in range( nsegments ): segments[i] *= conv[i]
        return segments
    
    def calculate_sfactor_rwa( self, params, proc = 0 ):
        nsegments = int(self.clients[proc].communicate( "UPDATE_SEGMENTS_RWA", params ))
        segments = [ self.clients[proc].communicate( "GET_CALCULATED_SEGMENT", [i] ) for i in range( nsegments ) ]
        conv = [ self.clients[proc].communicate( "GET_CALCULATED_CONV", [i] ) for i in range( nsegments ) ]
        for i in range( nsegments ): segments[i] *= conv[i]
        return segments
    
    def extrap_mode( self ):
        def _set_extrap(client):
            client.communicate( "SET_EXTRAP_MODE", [0] )
            client.communicate( "INITIALIZE", [0] )
        with ThreadPoolExecutor(max_workers=len(self.clients)) as ex:
            ex.map(_set_extrap, self.clients)

    def data_mode( self ):
        def _set_data(client):
            client.communicate( "SET_DATA_MODE", [0] )
            client.communicate( "INITIALIZE", [0] )
        with ThreadPoolExecutor(max_workers=len(self.clients)) as ex:
            ex.map(_set_data, self.clients)