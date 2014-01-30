import subprocess 
class NotEnergy(Exception):
    pass

class MakeGetEnergyFunc:
    """class for output file parsing, if en(str) is appropriate function, that returns total energy from string, then
    the MakeGetEnergyGunc(en) is function that returns Total energy for the given output Stream"""
    def __init__(self,get_energy):
        """get_energy -- callable that takes string and returns total energy if string contains it and nothing otherwise"""
        self.get_energy=get_energy
    def __call__(self,out_file):
        with open(out_file,'r') as out_stream:
            for line in out_stream:
                res =  self.get_energy(line)
                if res is not None: 
                    return res
            #Error 
        if res is None:
           raise NotEnergy

class MakeRunProgram:
    """class for running computation program recieves string with command, that
    executes program as prog_str and function that generates input file for the
    program as transform_input.  
    
    The transform_input function has to has following
    signature: res = transform_input(*args), where type of res is the
    string.

    MakeRunProgram is the callable object that passes it args to the transform_input function and
    runs shell command defined by the string run_str concatenated from the prog_str and res strings.

    Usage:
    Suppose we want to make function that runs the program "prg" with input
    file "inp.inp", we would run that program manually in the shell as follows:
    $>prg inp.inp.

    we have to define function "transform_input_for_prg" with following form:
    
    def transform_input_for_prg(arg1,arg2,...,argn):
        ... #writing to the inp.inp something depends on arg1,...,argn
        return "inp.inp"

    Then the function that we wanted is
    fun = MakeRunProgram("prg",transform_input_for_prg)

    """
    def __init__(self, prog_str, transform_input):
        self.prog_str = prog_str #the name of command to execute
        self.transform_input = transform_input #takes input file_name transforms it with accounting for args and returns command suffix
    def __call__(self,*args):
        run_str = self.prog_str + " " + self.transform_input(*args) 
        print ('"{}"'.format(run_str))
        subprocess.call(run_str,shell = True)

def dirac_en(str):
    if 'Total energy' in str:
        tmp = str.split(':')
        return float(tmp[1])

def transform_input_diatomic(in_file,true_in_file,dist):
    import re
    with open(in_file,'r') as in_stream:
        with open(true_in_file,'w') as out_stream:
            for line in in_stream:
               out_stream.write(re.sub(r'!dist!',str(dist),line)) #stream must contain !dist! in place of some coordinate (z for example) of some atom
    return '--mol="{}"'.format(true_in_file)

def transform_input_eval(in_file,true_in_file,d):
""" transforms input in the following way:
    All strings remains unchanged except those, that contain !...! substring. The ... is some python code. In this
    code all builtin python functions as well  math module function sqrt and argument d can be accessed. The substring !...! is replaced with result of this code. 
    Example:
    The string !d*0.1! with value of the d argument equal to float number 2 will be replaced to the 0.2
"""
    import re
    import math
    def replace_func(matchObj,d=d): #d=d is magic for eval and re.sub combination
        return str(eval(matchObj.group(1),globals(),{'sqrt':math.sqrt,'d':d}))
    evalp = re.compile(r'!([^!]+)!') #eval_pattern
    with open(in_file,'r') as in_stream:
        with open(true_in_file,'w') as out_stream:
            for line in in_stream:
               out_line = re.sub(evalp,replace_func,line)
               out_stream.write(out_line)
    return '--mol="{}"'.format(true_in_file)
