using JuMP, Gurobi, Random

function formulationMTZ(relax::Bool=false)
  include("instance.txt")
  println("MTZ******************************* n=", n, " m=", m, " p=", p)
  #--------------------------------------------------------------
  #declaration du modele
model = Model(Gurobi.Optimizer)
model=Model(
        optimizer_with_attributes(
          Gurobi.Optimizer, "Presolve" => -1
        ))
        
set_optimizer_attribute(model, "TimeLimit", 3600)

#declaration des variables
@variable(model, 0 <= c[1:n] <=1)  # vaut 1 ssi client i couvert
@variable(model, t[j in 1:m] >= 0)  
if (relax==false) 
  @variable(model, 0 <= y[1:m], Bin)  #vaut 1 si un relais est place en position j
  @variable(model, 0 <= x[j in 1:m, j1 in 1:m ; j1!=j] , Bin)  # vaut 1 si l'arc (j,j1) fait partie de l'arborescence
else
  @variable(model, 0 <= y[1:m] <=1) 
  @variable(model, 0 <= x[j in 1:m, j1 in 1:m ; j1!=j] <=1) 
end

#declaration de l objectif
@objective(model, Max, sum(c[i]  for i in 1:n) )

#declaration des contraintes
@constraint(model, [i in 1:n], c[i] <=  sum(y[j] for j in 1:m if (d[i,j]<=Rcouv))  )
@constraint(model,  sum(y[j] for j in 1:m) == p)

@constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j && d_com[j,j1] > Rcom] , x[j,j1]== 0)
@constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j] ,  x[j,j1]+x[j1,j]<= y[j])
@constraint(model,  [j in 1:m] , sum( x[j1,j] for j1 in 1:m if j1!=j) <= y[j])
@constraint(model,  sum(x[j,j1] for j in 1:m for j1 in 1:m if j1!=j) == sum(y[j] for j in 1:m) -1)
@constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j] , t[j1] >= t[j] + 1 - p*(1-x[j,j1]))

#print(model
#resolution
JuMP.optimize!(model)   
#affichage des resultats
obj_value = JuMP.objective_value(model)
println("Objective value: ", obj_value)
println("positions des relais")
for j in 1:m
  if (JuMP.value(y[j])>0.5) 
    println(j)
  end
end
println("Arcs")
for j in 1:m for  j1 in 1:m    
    if (j1 !=j && JuMP.value(x[j,j1])>0.5) 
    println("j=", j, "j1=",  j1)
  end
end
end
end


################################################################################
################################################################################
function formulationMTZ2(relax::Bool=false)
include("instance.txt")
println("MTZ2******************************* n=", n, " m=", m, " p=", p)
#--------------------------------------------------------------
#declaration du modele
model = Model(Gurobi.Optimizer)
model=Model(
        optimizer_with_attributes(
          Gurobi.Optimizer, "Presolve" => -1
        ))
        
set_optimizer_attribute(model, "TimeLimit", 3600)

#declaration des variables
@variable(model, 0 <= c[1:n] <=1)  # vaut 1 ssi client i couvert
@variable(model, t[j in 0:m] >= 0)  
if (relax==false) 
  @variable(model, 0 <= y[1:m], Bin)  #vaut 1 si un relais est place en position j
  @variable(model, 0 <= x[j in 0:m, j1 in 1:m ; j1!=j] , Bin)  # vaut 1 si l'arc (j,j1) fait partie de l'arborescence
else
  @variable(model, 0 <= y[1:m] <=1) 
  @variable(model, 0 <= x[j in 0:m, j1 in 1:m ; j1!=j] <=1) 
end

#declaration de l objectif
@objective(model, Max, sum(c[i]  for i in 1:n) )

#declaration des contraintes
@constraint(model, [i in 1:n], c[i] <=  sum(y[j] for j in 1:m if (d[i,j]<=Rcouv))  )
@constraint(model,  sum(y[j] for j in 1:m) == p)
@constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j && d_com[j,j1] > Rcom] , x[j,j1]== 0)
@constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j] ,  x[j,j1]+x[j1,j]<= y[j])
@constraint(model,  [j in 1:m] , sum( x[j1,j] for j1 in 0:m if j1!=j) == y[j])
@constraint(model,  sum(x[j,j1] for j in 1:m for j1 in 1:m if j1!=j) == sum(y[j] for j in 1:m) -1)
@constraint(model,  sum(x[0,j] for j in 1:m ) == 1)
@constraint(model, [j in 1:m] , x[0,j] <= y[j])
#connexite
@constraint(model,  t[0] == 0)
@constraint(model,  [j in 0:m, j1 in 1:m ; j1!=j] , t[j1] >= t[j] + 1 - p*(1-x[j,j1]))

#print(model)

#resolution
JuMP.optimize!(model)  
#--------------------------------------------------------------

#affichage des resultats
obj_value = JuMP.objective_value(model)
println("Objective value: ", obj_value)
println("positions des relais")
for j in 1:m
  if (JuMP.value(y[j])>0.5) 
    println(j)
  end
end
println("Arcs")
for j in 1:m for  j1 in 1:m    
    if (j1 !=j && JuMP.value(x[j,j1])>0.5) 
    println("j=", j, "j1=",  j1)
  end
end
end
end

##################################################
function formulationFlot(relax::Bool=false)

  include("instance.txt")
  println("Flot******************************* n=", n, " m=", m, " p=", p)  
  #--------------------------------------------------------------
  #declaration du modele
  model = Model(Gurobi.Optimizer)
  model=Model(
          optimizer_with_attributes(
            Gurobi.Optimizer, "Presolve" => -1
          ))
          
  set_optimizer_attribute(model, "TimeLimit", 3600)
  #declaration des variables
  @variable(model, 0 <= c[1:n] <=1)  # vaut 1 ssi client i couvert
  if (relax==false) 
    @variable(model, 0 <= y[1:m], Bin)  #vaut 1 si un relais est place en position j
    @variable(model, 0 <= x[j in 0:m, j1 in 1:m ; j1!=j] , Bin)  # vaut 1 si l'arc (j,j1) fait partie de l'arborescence
  else
    @variable(model, 0 <= y[1:m] <=1) 
    @variable(model, 0 <= x[j in 0:m, j1 in 1:m ; j1!=j] <=1) 
  end
  @variable(model, g[j in 0:m, j1 in 1:m ; j1!=j] >= 0)  
  
  #declaration de l objectif
  @objective(model, Max, sum(c[i]  for i in 1:n) )
  
  #declaration des contraintes
  @constraint(model, [i in 1:n], c[i] <=  sum(y[j] for j in 1:m if (d[i,j]<=Rcouv))  )
  @constraint(model,  sum(y[j] for j in 1:m) == p)
  
  @constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j && d_com[j,j1] > Rcom] , x[j,j1]== 0)
  @constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j] ,  x[j,j1]+x[j1,j]<= y[j])
  @constraint(model,  [j in 1:m] , sum( x[j1,j] for j1 in 0:m if j1!=j) == y[j])
  @constraint(model,  sum(x[j,j1] for j in 1:m for j1 in 1:m if j1!=j) == sum(y[j] for j in 1:m) -1)
  @constraint(model,  sum(x[0,j1] for j1 in 1:m ) == 1)
  @constraint(model, [j in 1:m] , x[0,j] <= y[j])
  #connexite
  @constraint(model,  sum(g[0,j] for j in 1:m ) == p)
  @constraint(model,  [j in 1:m] , sum( g[j1,j] for j1 in 0:m if j1!=j)  - sum( g[j,j1] for j1 in 1:m if j1!=j)== y[j])
  @constraint(model,  [j in 0:m, j1 in 1:m ; j1!=j] ,  g[j,j1] <= p*x[j,j1])
  
  #print(model)
  #resolution
  JuMP.optimize!(model)
  #--------------------------------------------------------------
  #affichage des resultats
  obj_value = JuMP.objective_value(model)
  println("Objective value: ", obj_value)
  println("positions des relais")
  for j in 1:m
    if (JuMP.value(y[j])>0.5) 
      println(j)
    end
  end
  println("Arcs")
  for j in 1:m for  j1 in 1:m    
      if (j1 !=j && JuMP.value(x[j,j1])>0.5) 
      println("j=", j, "j1=",  j1)
    end
  end
  end
end

##################################################
function formulationMultiflot(relax::Bool=false)
  include("instance.txt")
  println("MultiFlot******************************* n=", n, " m=", m, " p=", p)  
#--------------------------------------------------------------
#declaration du modele
model = Model(Gurobi.Optimizer)
model=Model(
        optimizer_with_attributes(
          Gurobi.Optimizer, "Presolve" => -1
        ))
        
set_optimizer_attribute(model, "TimeLimit", 3600)

#declaration des variables
@variable(model, 0 <= c[1:n] <=1)  # vaut 1 ssi client i couvert
@variable(model, 0 <=f[j in 0:m, j1 in 1:m , k in 1:m; j1!=j]  <=1)  # =1 si (j,j1) est sur le chemin de 0 à k 

if (relax==false) 
  @variable(model, 0 <= y[1:m], Bin)  #vaut 1 si un relais est place en position j
  @variable(model, 0 <= x[j in 0:m, j1 in 1:m ; j1!=j] , Bin)  # vaut 1 si l'arc (j,j1) fait partie de l'arborescence
else
  @variable(model, 0 <= y[1:m] <=1) 
  @variable(model, 0 <= x[j in 0:m, j1 in 1:m ; j1!=j] <=1) 
end

#declaration de l objectif
@objective(model, Max, sum(c[i]  for i in 1:n) )

#declaration des contraintes
@constraint(model, [i in 1:n], c[i] <=  sum(y[j] for j in 1:m if (d[i,j]<=Rcouv))  )
@constraint(model,  sum(y[j] for j in 1:m) == p)

@constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j && d_com[j,j1] > Rcom] , x[j,j1]== 0)
@constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j] ,  x[j,j1]+x[j1,j]<= y[j])
@constraint(model,  [j in 1:m] , sum( x[j1,j] for j1 in 0:m if j1!=j) == y[j])
@constraint(model,  sum(x[j,j1] for j in 1:m for j1 in 1:m if j1!=j) == sum(y[j] for j in 1:m) -1)
@constraint(model,  sum(x[0,j] for j in 1:m ) == 1)
@constraint(model, [j in 1:m] , x[0,j] <= y[j])
#connexite
@constraint(model,  [k in 1:m] ,sum(f[0,j,k] for j in 1:m ) == y[k])
@constraint(model,  [k in 1:m] ,sum(f[j,k,k] for j in 0:m if j !=k)  - sum(f[k,j,k] for j in 1:m if j !=k)== y[k])

@constraint(model, [k in 1:m, j in 1:m; j!=k] , sum( f[j1,j, k] for j1 in 0:m if j1!=j)  - sum( f[j,j1, k] for j1 in 1:m if j1!=j)== 0)

@constraint(model,  [k in 1:m, j in 0:m, j1 in 1:m ; j1!=j] ,  f[j,j1, k] <= x[j,j1])
#@constraint(model,  [j in 1:m] ,  sum(f[0,j,k] for k in 1:m) == p*x[0,j])
#@constraint(model,  [j in 1:m] ,  f[0,j,j]  == x[0,j])
@constraint(model,  [ j in 0:m, j1 in 1:m ; j1!=j] ,  sum(f[j,j1, k] for k in 1:m) <= p*x[j,j1])
#print(model)

#resolution
JuMP.optimize!(model)
#--------------------------------------------------------------
#affichage des resultats
obj_value = JuMP.objective_value(model)
println("Objective value: ", obj_value)
println("positions des relais")
for j in 1:m
  if (JuMP.value(y[j])>0.5) 
    println(j)
  end
end
println("Arcs")
for j in 1:m for  j1 in 1:m    
    if (j1 !=j && JuMP.value(x[j,j1])>0.5) 
    println("j=", j, "j1=",  j1)
  end
end
end
end

##########################################################################
##########################################################################
function sepGenExpoModel(mode::String="sep_relax",maxCuts::Int64=50)
  include("instance.txt")
  println("Subtour******************************* n=", n, " m=", m, " p=", p)  
#--------------------------------------------------------------
  #declaration du modele
  model = Model(Gurobi.Optimizer)
  model=Model(
          optimizer_with_attributes(
            Gurobi.Optimizer, "Presolve" => -1
          ))
  #declaration des variables
  @variable(model, 0 <= z[1:n] <=1)  # vaut 1 ssi point i couvert
  @variable(model, 0 <= y[1:m] <=1)  #vaut 1 si est relais est place en position j
  @variable(model, 0<= x[j in 1:m, j1 in 1:m ; j1!=j]<=1)  #les relais j1 et j communiquent directement
  #@variable(model, t[j in 1:m] >= 0)  
  
  #declaration de l objectif
  @objective(model, Max, sum(z[i]  for i in 1:n) )
  
  #declaration des contraintes
  @constraint(model, [i in 1:n], z[i] <=  sum(y[j] for j in 1:m if (d[i,j]<=Rcouv))  )
  @constraint(model,  sum(y[j] for j in 1:m) == p)
  
  @constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j && d_com[j,j1] > Rcom] , x[j,j1]== 0)
  @constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j] , x[j,j1]+x[j1,j]<= y[j])
  @constraint(model,  [j in 1:m] , sum( x[j1,j] for j1 in 1:m if j1!=j) <= y[j])
  @constraint(model,  sum(x[j,j1] for j in 1:m for j1 in 1:m if j1!=j) == sum(y[j] for j in 1:m) -1)
  
  #print(model)

  #Solve model - with time limit and separate subtour constraints
    set_optimizer_attribute(model, "TimeLimit", 3600)
    set_silent(model)
  
    isOptimal = false
    best_lower_bound=0
    nbCuts = 0
    resTime = time()
  
    while (!isOptimal) && (nbCuts < maxCuts)
      optimize!(model)     
      isOptimal = true
   
      #Check existence of solution - check whether everything was solved properly
      if (termination_status(model) == MOI.OPTIMAL)
        println("Solved to optimality. Value : ",objective_value(model))
        best_lower_bound=objective_value(model)
      elseif (termination_status(model) == MOI.TIME_LIMIT && has_values(model))
        println("Time limit reached, primal solution available. Value : ",objective_value(model))
      elseif (termination_status(model) == MOI.INFEASIBLE)
        error("Problem infeasible, verify constraints or instance before solving.")
      else
        error("Problem could not be solved.")
      end
  
      #Separate subtour elimination constraints
      oldstd = stdout
      redirect_stdout(open("fich.log", "w"))
      exist_violated_cut= false
      
        sep = Model(Gurobi.Optimizer)
        redirect_stdout(oldstd)
        set_silent(sep)
        @variable(sep, a[1:m], Bin)
        @variable(sep, h[1:m], Bin)
        @variable(sep, w[j in 1:m, j1 in 1:m; j<j1] >= 0)
        @objective(sep, Max,  sum((value(x[j,j1])+value(x[j1,j]))*w[j,j1] for j in 1:m-1 for j1 in j+1:m)
               - sum(value(y[j])*a[j] for j in 1:m) + sum(value(y[j])*h[j] for j in 1:m)  )
        @constraint(sep, inf1[j=1:m-1,j1=j+1:m], w[j,j1] <= a[j])
        @constraint(sep, inf2[j=1:m-1,j1=j+1:m], w[j,j1] <= a[j1])
        @constraint(sep, inf3[j=1:m], h[j] <= a[j])
        @constraint(sep, sum(h[j] for j in 1:m) ==1)
        optimize!(sep)
        #println(objective_value(sep))
  
        #Check if some subtour constraint is not verified
        if (objective_value(sep) > 0.01)
          exist_violated_cut = true
          println("SEC violated. Adding constraint to master problem...")
          @constraint(model, sum((x[j,j1]+x[j1,j])*value(a[j])*value(a[j1]) for j in 1:m-1 for j1 in j+1:m)
            <= sum(value(a[j])*y[j] for j in 1:m) - sum(value(h[j])*y[j] for j in 1:m) )
          nbCuts += 1
        end
      
      isOptimal = !exist_violated_cut
      if isOptimal
        println("No SECs violated. Valid solution reached.")
      end
    end
    optimize!(model) 
    println("="^76)
    if (nbCuts >= maxCuts)
      println("Maximum number of iterations reached.")
    end
  
  
    #JuMP.value.(z)
    for j in 1:m ,j1 in 1:m    
      if (j1 !=j && JuMP.value(x[j,j1])>0) 
      println("j=", j, "j1=",  j1, "  ",JuMP.value(x[j,j1]) )
    end
  end
  
    println("----best lower bound = ", best_lower_bound)
    println("----number of added cuts = ", nbCuts)
    println("----Time spent up here : ",time()-resTime," s")
  
    if (true)
      for j in 1:m
        #set_integer(y[j])
        for j1 in 1:m if (j !=j1)
          set_integer(x[j,j1])
        end
        end
      end
      println("Re-optimizing with binary variables and MTZ for optimal solution...")
  
      @variable(model, t[j in 1:m] >= 0)  
      @constraint(model,  [j in 1:m, j1 in 1:m ; j1!=j] , t[j1] >= t[j] + 1 - p*(1-x[j,j1]))
  
      unset_silent(model)
      optimize!(model)
      println("fin opt")
      #Check existence of solution - check whether everything was solved properly
      if (termination_status(model) == MOI.OPTIMAL)
          println("Solved to optimality. Value : ",objective_value(model))
        elseif (termination_status(model) == MOI.TIME_LIMIT && has_values(model))
          println("Time limit reached, primal solution available. Value : ",objective_value(model))
        elseif (termination_status(model) == MOI.INFEASIBLE)
          error("Problem infeasible, verify constraints or instance before solving.")
        else
          error("Problem could not be solved.")
        end
      end
  
    #println("Best value found : ",objective_value(model))
    println("Total time spent : ",time()-resTime," s")
  # -----------------------------------------------------------
  #affichage des resultats
  obj_value = JuMP.objective_value(model)
  println("Objective value: ", obj_value)
  println("positions des relais")
  
  for j in 1:m
    if (JuMP.value(y[j])>0.5) 
      println(j)
    end
  end
  println("Arcs")
  for j in 1:m for  j1 in 1:m    
      if (j1 !=j && JuMP.value(x[j,j1])>0.5) 
      println("j=", j, "j1=",  j1)
    end
  end
  end
  end
  