import math
def aire(tri):
    #assert len(tri)==3, "trois cotes svp"
    if(len(tri))!=3:
        print ("trois cotes svp")
        return None
    assert tri[0]>0 and tri[1]>0 and tri[2]>0, "cotes positifs svp"
    #assert for i in tri:
     #          tri[:]-2*i >0
      #      , "somme de deux cotes toujours plus grande que le dernier cote"
    """calcul de l'aire d'un triangle de cotes abc"""
    a,b,c=tri
    cos_al=(b**2+c**2-a**2)/(2*b*c)
    al=math.acos(cos_al)
    h=b*math.sin(al)
    return h*c/2

