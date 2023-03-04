import sys                                                                       
import qxelarator                                                                       
qx = qxelarator.QX()                                                                       
qx.set(sys.argv[1])                                                                       
json_string = qx.execute()                                                                       
print(json_string)      