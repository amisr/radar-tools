#! /usr/bin/env python

"""
xxxxx

~M. Nicolls
last revised: xx/xx/2007

~P. Reyes
changes:
 - in read_whole_h5file include an exclusion list use open_file in tables.
 - in read_partial_h5file there is a groups to do and groups not to do list entry.
 - the data numpy arrays are stored as a inherited array class that allows
   for metadata containing the attributes.
 - 22Nov2017 by P. Reyes:
   fixed the groupstodo groupsnottodo behaviour in print_info_dict
 - 27Nov2017 by P. Reyes:
   in read_partial_h5file it is possible to just read an element of a group
   e.g. /Setup/BeamcodeMap will only load BeamcodeMap from /Setup group
 - 11Dec2017 by P. Reyes:
   In read_whole_h5file groupsnottodo items need to include also '/'
 - 18Dec2017 by P. Reyes:
   In print_info_dict, also sub-parts of primary keys are tested for groupsnottodo
   In print_info_dict, add print_array parameter to print the array contents. 
   Also avoid printing the configuration files, but include the size of them.
"""
from __future__ import print_function
import tables
import os
import numpy as np

class __dummy__:
    pass
class DotAccess:
    """
    Accessing the data using dot.<TAB> using nested  classes.
    """
    def __init__(self,file0):
        self.contents_dict = file0
        for key,val in file0.items():
            keys = key.split('/')
            keys[0] = 'self'
            for key0 in keys:
                key0 = key0.replace('.','_dot_') # replace names that have '.'
            for i,group in enumerate(keys[1:]):
                if len(group)>0:
                    self.__add__('.'.join(keys[:i+1]),group)

            for key1,val2 in val.items():
                key1 = key1.replace('.','_dot_') # replace names that have '.'
                if len(keys[-1])>0:
                    exec("{0}.{1}=val2".format('.'.join(keys),key1))
                else:
                    exec("{0}.{1}=val2".format('.'.join(keys[:-1]),key1))
    def __add__(self,class0,element):
        if not eval("hasattr({0},'{1}')".format(class0,element)):
            toexec = "{0}.{1} = __dummy__()".format(class0,element)
            exec(toexec)

    def print_info(self,groupstodo=[],groupsnottodo=[],print_array=False,plotattrs=False,
                    avoid_attrs=['CLASS','TITLE','VERSION']):
        print_info_dict(self.contents_dict,groupstodo=groupstodo,
                groupsnottodo=groupsnottodo,print_array=print_array,
                plotattrs=plotattrs,avoid_attrs=avoid_attrs)

class AArray(np.ndarray):
    """Array with attributes."""

    def __new__(cls, array, attrs):
        obj = np.asarray(array).view(cls)
        obj.attrs = attrs
        obj.attributes = obj.attrs
        return obj

    def __array_finalize__(self, obj):
        if obj is None: return
        self.attrs = getattr(obj, 'attr', None)

def read_whole_h5file(fname,groupsnottodo=[]):

    h5file=tables.open_file(fname)
    output={}
    for group in h5file.walk_groups("/"):
        if group._v_pathname in groupsnottodo:
            print("skipping",group._v_pathname)
            continue
        output[group._v_pathname]={}
        for array in h5file.list_nodes(group, classname = 'Array'):
            attrs = {}
            for attr in [x for x in array.attrs.__dict__.keys() if x[0]!='_']:
                # getting the attribute names that don't start with underscore
                attrs.update({attr:array.attrs.__dict__[attr]})
            output[group._v_pathname][array.name]=AArray(array.read(),attrs)
        for table in h5file.list_nodes(group, classname = 'Table'):
            attrs = {}
            for attr in [x for x in table.attrs.__dict__.keys() if x[0]!='_']:
                # getting the attribute names that don't start with underscore
                attrs.update({attr:table.attrs.__dict__[attr]})
            output[group._v_pathname][table.name]=AArray(table.read(),attrs)
#            try:
#                print array.Size
#            except:
#                pass
    h5file.close()

    return DotAccess(output)

def read_partial_h5file(fname,groupstodo=[],groupsnottodo=[]):

    h5file=tables.open_file(fname)
    output={}
    for group in h5file.walk_groups("/"):
        bout = True
        partialgroup = False
        if len(groupstodo)>0:
            #print ".......checking",group._v_pathname
            if group._v_pathname not in groupstodo:
                # try to find the case of a sub group
                start = 0
                start = group._v_pathname.find('/',start+1)
                bout = False
                while start>0:
                    if group._v_pathname[0:group._v_pathname.find('/',start)] in groupstodo:
                        bout = True
                        break
                    start = group._v_pathname.find('/',start+1) # search for next level of '/'
                if not bout:
                    # means that the group is not part of a request
                    # Do a search of base group matching:

                    for grouptodo in groupstodo:
                        start = 0
                        start = grouptodo.find('/',start+1)
                        while start>0:
                            if grouptodo[0:grouptodo.find('/',start)] == group._v_pathname:
                                # if the group is part of any of the groupstodo list
                                bout = True
                                partialgroup = True
                                break
                            start = grouptodo.find('/',start+1) # search for next level of '/'
        if group._v_pathname in groupsnottodo:
            bout = False
        else:
            # try to find the case of a sub group
            start = 0
            start = group._v_pathname.find('/',start+1)
            while start>0:
                if group._v_pathname[0:group._v_pathname.find('/',start)] in groupsnottodo:
                    bout = False
                start = group._v_pathname.find('/',start+1)
        if bout:
            output[group._v_pathname]={}
            for array in h5file.list_nodes(group, classname = 'Array'):
                attrs = {}
                for attr in [x for x in array.attrs.__dict__.keys() if x[0]!='_']:
                    # getting the attribute names that don't start with underscore
                    attrs.update({attr:array.attrs.__dict__[attr]})
                if partialgroup:
                    if group._v_pathname + '/'+ array.name in groupstodo:
                        output[group._v_pathname][array.name]=AArray(array.read(),attrs)
                else:
                    output[group._v_pathname][array.name]=AArray(array.read(),attrs)
    h5file.close()

    return DotAccess(output)


def print_info_dict(indata,groupstodo=[],groupsnottodo=[],print_array=False,plotattrs=False,
        avoid_attrs=['CLASS','TITLE','VERSION']):
    """
    avoid_attrs: list of attributes to skip printing
    """
    if not plotattrs:
        print("Default argument: plotattrs=False")
    import numpy as np
    def pt(vtype,nlen):
        return vtype+(nlen-len(vtype))*'.'
    def print_line(key1,vdtype,typelen,vshape,val):
        if val.nbytes > (1024.)**2:
            memstr = "%8.3g MB"%(val.nbytes/(1024.)**2)
        elif val.nbytes > 1024.:
            memstr = "%8.3g kB"%(val.nbytes/1024.)
        else:
            memstr = "%8.3g Bytes"%(val.nbytes)

        print_1val = False
        if val.nbytes <= 16:
            print_1val = True
        if val.dtype.type == np.string_ :
            if val.size == 1:
                if val.item().find(b'\n')<0:
                    print_1val = True
        if print_1val:
            print("  |__",key1,(18-len(key1))*'.',pt(vdtype,typelen),pt(vshape,20),memstr,val)
        else:
            val_squeeze = val.squeeze()
            if val_squeeze.shape.__len__() == 1:
                try:
                    if np.iscomplex(val_squeeze[0]):
                        first_val = "[0]=(%g+%gj)"%(val_squeeze[0].real,
                                val_squeeze[0].imag)
                        last_val = "[-1]=(%g+%gj)"%(val_squeeze[-1].real,
                                val_squeeze[-1].imag)
                    else:
                        first_val = "[0]=%g"%val_squeeze[0]
                        last_val = "[-1]=%g"%val_squeeze[-1]
                    print("  |__",key1,(18-len(key1))*'.',pt(vdtype,typelen),pt(vshape,20),memstr,\
                        first_val,last_val)
                except:
                    print("  |__",key1,(18-len(key1))*'.',pt(vdtype,typelen),pt(vshape,20),memstr)
            else:
                print("  |__",key1,(18-len(key1))*'.',pt(vdtype,typelen),pt(vshape,20),memstr)
            if print_array:
                print(val)
    def print_str(key1,vtypename,typelen,val):
        if len(val) > 1024**2:
            memstr = "%8.3g MB"%(len(val)/(1024.)**2)
        if len(val) > 1024:
            memstr = "%8.3g kB"%(len(val)/1024.)
        else:
            memstr = "%8.3g Bytes"%(len(val))
        print("  |__",key1,(18-len(key1))*'.',vtypename,(typelen-len(vtypename))*'.',memstr)
    typelen = 13
    for key0 in sorted(indata.keys()):
        if len(groupstodo)>0:
            bprint = False
            for grouptodo in groupstodo:
                if grouptodo == key0:
                    bprint = True
                    break
                elif grouptodo[0:grouptodo.find('/',1)] == key0:
                    bprint = True
                    break
            if not bprint:
                continue
        if len(groupsnottodo)>0:
            bprint = True
            for i in range(1,len(key0.split('/'))):
                "searching also subgroups"
                if "/".join(key0.split(('/'))[:i+1]) in groupsnottodo:
                    bprint = False
                    break
            if not bprint:
                continue
        print(key0)
        for key1 in sorted(indata[key0].keys()):
            if len(groupsnottodo)>0:
                if key0+'/'+key1 in groupsnottodo:
                    continue
            if len(groupstodo)>0:
                if key0 not in groupstodo:
                    if key0+'/'+key1 not in groupstodo:
                        continue
            val = indata[key0][key1]
            vtype = type(val)
            vtypename = type(val).__name__
            if vtype in [np.ndarray,AArray]:
                vshape = val.shape.__str__()
                vdtype = val.dtype.name
                #if val.attrs['CLASS']=="TABLE":
                #    print "  |__",key1,(18-len(key1))*'.'
                #else:
                print_line(key1,vdtype,typelen,vshape,val)
            elif vtype in [float,int]:
                print("  |__",key1,(18-len(key1))*'.',vtypename,(typelen-len(vtypename))*'.',val)
            elif vtype in [str]:
                if "\n" in val:
                    print_str(key1,vtypename,typelen,val)
                    start = 0
                    for i in range(5):
                        end = val.find('\n',start)
                        print(8*" ",val[start:end])
                        start=end+1
                    print(8*" ","...")
                else:
                    print("  |__",key1,(18-len(key1))*'.',vtypename,(typelen-len(vtypename))*'.',val)
            else:
                print("  |__",key1,(18-len(key1))*'.',vtypename,(typelen-len(vtypename))*'.')
            if vtype == AArray and plotattrs:
                for key2 in sorted(val.attrs.keys()):
                    if key2 in avoid_attrs:
                        continue
                    if type(plotattrs)==list:
                        if key2 not in plotattrs:
                            continue
                    print(4*" "+"|_",key2,(11-len(key2))*'.',val.attrs[key2])



def write_outputfile(fhandle,dict2do,keys2do=[],groupname='',name='',grouploc='/'):
    # Tested with pytables 2.0.1

    if groupname == '':
        group=fhandle.root
    else:
        if fhandle.__contains__(grouploc+groupname):
            group='/'+groupname
        else:
            #fhandle.root
            group=fhandle.create_group(grouploc, groupname, 'Dataset')

    if len(keys2do)==0:
        try:
            fhandle.remove_node(group,name)
        except:
            ''
        fhandle.create_array(group,name, dict2do, "Dataset")
    else:
        for key in keys2do:
            if type(dict2do[key]) is dict:
                write_outputfile(fhandle,dict2do[key],keys2do=dict2do[key].keys(),groupname=key,grouploc='/'+group._v_name)
            else:
                fhandle.create_array(group, key, dict2do[key], "Dataset")

    return

def createh5groups(fhandle,h5Groups):
    # creates groups
    for group in h5Groups:
        gp,gn = os.path.split(group[0])
        fhandle.createGroup(gp,gn,group[1])
    return

def createStaticArray(fhandle,path,data,keys2do=[]):
    # creates a static array
    if len(keys2do)==0:
        dp,dn = os.path.split(path)
        fhandle.createArray(dp,dn,data,'Static array')
    else:
        for key in keys2do:
            fhandle.createArray(path,key,data[key],'Static array')
    return

def setAtrributes(fhandle,data):
    for key in data.keys():
        for attr in data[key]:
            try:  fhandle.setNodeAttr(key,attr[0],attr[1])
            except: ''
    return

