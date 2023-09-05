# -*- coding: utf-8 -*-
import os
import sys
import logging
import numpy as np
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

try:
    from openiam import SystemModel, ComponentModel
except ImportError as err:
    print('Unable to load IAM class module: {}'.format(err))


class Stratigraphy(ComponentModel):
    """
    The Stratigraphy component is a component containing parameters describing
    the structure of the storage reservoir system. The stratigraphy component
    allows the number of shale (or aquitard) layers to be set, thus, setting the
    total number of layers in the system. Each shale or aquifer layer can
    take on the default thickness for that layer type or be assigned a user defined value
    with ``shale#Thickness`` or ``aquifer#Thickness`` keywords where # is replaced
    by an index of the layer (e.g., ``shale1Thickness``). Layers are numbered from
    the bottom upward: shale 1 is the layer above the storage reservoir,
    and, with N shale layers total, shale N is the surface layer.

    .. image:: ../../images/ShaleLayers.png
       :align: center
       :scale: 100%
       :alt: Layer ordering

    The description of the component's parameters is provided below.

    * **numberOfShaleLayers** [-] (3 to 30) - number of shale layers in the
      system (default: 3). The shale units must be separated by an aquifer.

    * **shaleThickness** [|m|] (1 to 1600) - thickness of shale layers (default
      250). Thickness of shale layer 1, for example, can be defined by
      **shale1Thickness**; otherwise, shale layers for which the thickness is not
      defined will be assigned a default thickness.

    * **aquiferThickness** [|m|] (1 to 1600) - thickness of aquifers (default: 100).
      Thickness of aquifer 1, for example, can be defined by **aquifer1Thickness**;
      otherwise, aquifers for which the thickness is not defined will be assigned
      a default thickness.

    * **reservoirThickness** [|m|] (1 to 1600) - thickness of reservoir (default: 50)

    * **datumPressure** [|Pa|] (80,000 to 300,000) - pressure at the top of the
      system (default: 101,325)

    * **depth** [|m|] (5 to 30,000) - depth to the top of reservoir (default: 950).

    For control file examples the following composite parameters are produced:

    * **shale#Depth** [|m|] (boundaries depend on the user input) - depth
      to the base of the shale layer with index # (default value is not defined)

    * **aquifer#Depth** [|m|] (boundaries depend on the user input) - depth
      to the base of the aquifer layer with index # (default value is not defined)

    * **reservoirDepth** [|m|] (boundaries depend on the user input) - depth
      to the base of the reservoir (default value is not defined).

    For script examples these parameters have to be added explicitly as composite
    parameters of the stratigraphy component.
    """
    def __init__(self, name, parent):
        """
        Constructor method of Container Class.
        """
        super().__init__(name, parent, model='')

        # Add type attribute
        self.class_type = 'Stratigraphy'

        # Set default parameters of the component model
        self.add_default_par('numberOfShaleLayers', value=3)
        self.add_default_par('shaleThickness', value=250.0)
        self.add_default_par('aquiferThickness', value=100.0)
        self.add_default_par('reservoirThickness', value=50.0)
        self.add_default_par('datumPressure', value=101325.0)
        self.add_default_par('depth', value=950.0)   # depth to the top of the reservoir

        # Define dictionary of boundaries
        self.pars_bounds = dict()
        self.pars_bounds['numberOfShaleLayers'] = [3, 30]
        self.pars_bounds['shaleThickness'] = [1.0, 1600.0]
        self.pars_bounds['aquiferThickness'] = [1.0, 1600.0]
        self.pars_bounds['reservoirThickness'] = [1.0, 1600.0]
        self.pars_bounds['datumPressure'] = [8.0e+4, 3.0e+5]
        self.pars_bounds['depth'] = [5.0, 30000.0]

        # Indicate that the component should not be run
        self.default_run_frequency = 0
        self.run_frequency = 0

        debug_msg = 'Stratigraphy component created with name {}'.format(name)
        logging.debug(debug_msg)

    def check_input_parameters(self, p):
        """
        Check whether input parameters fall within specified boundaries.

        :param p: input parameters of component model
        :type p: dict
        """
        debug_msg = 'Input parameters of component {} are {}'.format(self.name, p)
        logging.debug(debug_msg)

        for key, val in p.items():
            if not key in self.gridded_pars:
                warn_msg = ''.join([
                    'Parameter {} of Stratigraphy component {} ',
                    'is out of boundaries.']).format(key, self.name)
                if (key[0:5] == 'shale' and key[-9:] == 'Thickness'):
                    if (val < self.pars_bounds['shaleThickness'][0]) or (
                            val > self.pars_bounds['shaleThickness'][1]):
                        logging.warning(warn_msg)
                elif key[0:7] == 'aquifer' and key[-9:] == 'Thickness':
                    if (val < self.pars_bounds['aquiferThickness'][0]) or (
                            val > self.pars_bounds['aquiferThickness'][1]):
                        logging.warning(warn_msg)
                elif key[0:5] == 'depth':
                    if (val < self.pars_bounds['depth'][0]) or (
                            val > self.pars_bounds['depth'][1]):
                        logging.warning(warn_msg)
                elif key in self.pars_bounds:
                    if (val < self.pars_bounds[key][0]) or (
                            val > self.pars_bounds[key][1]):
                        logging.warning(warn_msg)
                else:
                    warn_msg = ''.join([
                        'Parameter {} is not recognized as an input parameter ',
                        'of Stratigraphy component {}.']).format(key, self.name)
                    logging.warning(warn_msg)

    def connect_with_system(self):
        """
        Code to add stratigraphy to system model for control file interface.
        """
        # Check if numberOfShaleLayers is among stochastic parameters
        if 'numberOfShaleLayers' in self.pars:
            err_msg = ''.join(['Parameter numberOfShaleLayers cannot be ',
                               'stochastic for the control file interface.'])
            logging.error(err_msg)
            raise TypeError(err_msg)
        elif 'numberOfShaleLayers' in self.deterministic_pars:
            numberOfShaleLayers = self.deterministic_pars['numberOfShaleLayers'].value
        else:
            warn_msg = ''.join([
                'Parameter numberOfShaleLayers is not defined ',
                'in the control file interface. The parameter ',
                'will be assigned a default value of 3.'])
            logging.warn(warn_msg)
            numberOfShaleLayers = self.default_pars['numberOfShaleLayers'].value

        # Check whether all other needed parameters are defined by user
        par_names = ['shale{}Thickness'.format(ind) \
                     for ind in range(1, numberOfShaleLayers+1)] + [
                         'aquifer{}Thickness'.format(ind) \
                             for ind in range(1, numberOfShaleLayers)]

        for ind, par_nm in enumerate(par_names):
            if (par_nm not in self.pars) and (par_nm not in self.deterministic_pars):
                if ind < numberOfShaleLayers:
                    default_value = self.default_pars['shaleThickness'].value
                else:
                    default_value = self.default_pars['aquiferThickness'].value
                warn_msg = ''.join([
                    'Parameter {} is not defined in the control file ',
                    'interface. The parameter will be assigned ',
                    'a default value of {}.']).format(par_nm, default_value)
                logging.warn(warn_msg)
                self.add_par(par_nm, value=default_value, vary=False)

        # Add composite parameters of the stratigraphy component
        # representing depth to the bottom of each shale layer
        for ind in range(1, numberOfShaleLayers+1):
            par_expr = '+'.join(
                ['{}.shale{}Thickness'.format(self.name, j) for j in range(
                    ind, numberOfShaleLayers+1)]+[
                        '{}.aquifer{}Thickness'.format(self.name, j) for j in range(
                            ind, numberOfShaleLayers)])
            self.add_composite_par('shale{}Depth'.format(ind), par_expr)

        # Add composite parameters of the stratigraphy component
        # representing depth to the bottom of each aquifer layer
        for ind in range(1, numberOfShaleLayers):
            par_expr = '+'.join(
                ['{}.shale{}Thickness'.format(self.name, j) for j in range(
                    ind+1, numberOfShaleLayers+1)]+[
                        '{}.aquifer{}Thickness'.format(self.name, j) for j in range(
                            ind, numberOfShaleLayers)])
            self.add_composite_par('aquifer{}Depth'.format(ind), par_expr)

        # Add composite parameter of the stratigraphy component
        # representing depth to the bottom of reservoir
        par_expr = '+'.join(
            ['{}.shale{}Thickness'.format(self.name, j) for j in range(
                1, numberOfShaleLayers+1)]+[
                    '{}.aquifer{}Thickness'.format(self.name, j) for j in range(
                        1, numberOfShaleLayers)]+['{}.reservoirThickness'.format(self.name)])
        self.add_composite_par('reservoirDepth', par_expr)


if __name__ == "__main__":
    try:
        from openiam import SimpleReservoir
    except ImportError as err:
        print('Unable to load IAM class module: '+str(err))

    logging.basicConfig(level=logging.WARNING)
    # Define keyword arguments of the system model
    num_years = 50
    time_array = 365.25*np.arange(0.0, num_years+1) # time is in days
    sm_model_kwargs = {'time_array': time_array}
    sm = SystemModel(model_kwargs=sm_model_kwargs)

    strata = sm.add_component_model_object(Stratigraphy(name='strata', parent=sm))

    # Add parameters of reservoir component model
    strata.add_par('numberOfShaleLayers', value=3, vary=False)
    strata.add_par('shale1Thickness', min=30.0, max=50., value=40.0)
    strata.add_par('shale2Thickness', min=40.0, max=60., value=50.0)

    # Add reservoir component
    sres = sm.add_component_model_object(SimpleReservoir(name='sres', parent=sm))

    # Add parameters of reservoir component model
    sres.add_par_linked_to_par('numberOfShaleLayers',
                               strata.deterministic_pars['numberOfShaleLayers'])
    sres.add_par('injRate', min=0.4, max=0.6, value=0.5)
    sres.add_par_linked_to_par('shale1Thickness', strata.pars['shale1Thickness'])
    sres.add_par_linked_to_par('shale2Thickness', strata.pars['shale2Thickness'])
    sres.add_par_linked_to_par('shale3Thickness',
                               strata.default_pars['shaleThickness'])

    sres.add_par_linked_to_par('aquifer1Thickness',
                               strata.default_pars['aquiferThickness'])
    sres.add_par_linked_to_par('aquifer2Thickness',
                               strata.default_pars['aquiferThickness'])
    sres.add_par_linked_to_par('aquifer3Thickness',
                               strata.default_pars['aquiferThickness'])

    sres.add_par_linked_to_par('reservoirThickness',
                               strata.default_pars['reservoirThickness'])

    sres.add_par_linked_to_par('datumPressure',
                               strata.default_pars['datumPressure'])

    sres.add_obs('pressure')

    # Run system model using current values of its parameters
    sm.forward()

    print('------------------------------------------------------------------')
    print('                  Forward method illustration ')
    print('------------------------------------------------------------------')

    # Print pressure
    print('Pressure', sm.collect_observations_as_time_series(sres, 'pressure'),
          sep='\n')
