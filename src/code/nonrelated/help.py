print("Sniga S niga: ", end="")
rinda = int(input())

counter = 1  
rindus = rinda  

for i in range(rinda): 
    
    for k in range(rindus - 1): 
        print(" ", end="")  
    
  
    for j in range(counter): 
        print("*", end="")  
    
    print()  
    
    counter += 2  
    rindus -= 1 

    
