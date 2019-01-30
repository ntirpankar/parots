function result_process_PowerModel(result1)
    V=zeros(nbus);
    Pg=zeros(ngen);
    Qg=zeros(ngen);
    for i=1:nbus

        V[i]=result1["solution"]["bus"][string(i)]["vm"]
        #Theta[i]=result["solution"]["bus"][string(i)]["va"]
    end
    for i=1:ngen
        Pg[i]=result1["solution"]["gen"][string(i)]["pg"]
        Qg[i]=result1["solution"]["gen"][string(i)]["qg"]
        #Theta[i]=result["solution"]["bus"][string(i)]["va"]
    end
    #println(Pg)
    return(Pg, V)
end
function result_process_PowerModel_PF(result2)
    V=zeros(nbus);
    Theta=zeros(nbus);
    for i=1:nbus
        V[i]=result2["solution"]["bus"][string(i)]["vm"]
        Theta[i]=result2["solution"]["bus"][string(i)]["va"]
    end
    Pg=zeros(ngen)
    Qg=zeros(ngen)
    for i=1:ngen
        if (NaN in Set(result2["solution"]["gen"][string(i)]["pg"]))
            Pg[i]=0;
        else
            Pg[i]=result2["solution"]["gen"][string(i)]["pg"]
        end
        if (NaN in Set(result2["solution"]["gen"][string(i)]["qg"]))
            Qg[i]=0;
        else
            Qg[i]=result2["solution"]["gen"][string(i)]["qg"]
        end
        #Pg[i]=result2["solution"]["gen"][string(i)]["pg"]
        #Qg[i]=result2["solution"]["gen"][string(i)]["qg"]
        #Theta[i]=result["solution"]["bus"][string(i)]["va"]
    end
    #println(Pg)
    #return(V, Theta)
    return(Pg, Qg, V, Theta)
end
