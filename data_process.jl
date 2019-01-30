function data_process(network, t)
    global ncont=0;
    global SetCont=[1 2 3 4 5 6 8 10 11 12];
    #global network=PowerModels.parse_file(case)
    global nbus = length(network["bus"])
    global nbr = length(network["branch"])
    global ngen = length(network["gen"])
    global nshunt = length(network["shunt"])
    global linearSegments = 20;
    global linearCost=zeros(ngen,linearSegments)
    global lin_size=zeros(ngen)
    global cMin=zeros(ngen)
    global bus_array=zeros(nbus)

    #bus_array=map(x->(v = tryparse(Int64,x); isnull(v) ? 0.0 : get(v)), collect(keys(network["bus"])))
    for i=1:ngen
        #println(find(x->x==network["gen"][string(i)]["gen_bus"], bus_array)[1])
    end
     #network;
    #=
    if(t==1)
        network_temp=copy(network);
        network_temp["bus"]=Dict()
        global bus_array=Array(Int64)
        bus_array=map(x->(v = tryparse(Int64,x); isnull(v) ? 0.0 : get(v)), collect(keys(network["bus"])))
        for i=1:nbus
            network_temp["bus"][string(i)]=network["bus"][string(bus_array[i])]
        end
        #network_temp["gen"]=Dict()
        for i=1:ngen
            network_temp["gen"][string(i)]["gen_bus"]=find(x->x==network_temp["gen"][string(i)]["gen_bus"], bus_array)[1]
        end
        for i=1:nbr
            network_temp["branch"][string(i)]["f_bus"]=find(x->x==network_temp["branch"][string(i)]["f_bus"], bus_array)[1]
            network_temp["branch"][string(i)]["t_bus"]=find(x->x==network_temp["branch"][string(i)]["t_bus"], bus_array)[1]
        end
        for i=1:length(network["load"])
            network_temp["load"][string(i)]["load_bus"]=find(x->x==network_temp["load"][string(i)]["load_bus"], bus_array)[1]
        end
        for i=1: length(network["shunt"])
            network_temp["shunt"][string(i)]["shunt_bus"]=find(x->x==network_temp["shunt"][string(i)]["shunt_bus"], bus_array)[1]
        end
        network=copy(network_temp);
    end
    =#
    #=
    for i=1:ngen
        gen=network["gen"][string(i)]
        cMin[i]=gen["cost"][1]*gen["pmin"]^2+gen["cost"][2]*gen["pmin"]+gen["cost"][3];
        if (gen["pmax"]>gen["pmin"])
            lin_size[i]= (gen["pmax"]-gen["pmin"])/linearSegments;
            for l=1:linearSegments
                low=gen["pmin"]+lin_size[i]*(l-1);
                high=gen["pmin"]+lin_size[i]*(l);
                linearCost[i,l]=(gen["cost"][1]*high^2-gen["cost"][1]*low^2)/lin_size[i]+gen["cost"][2];
            end
        end
    end=#
    #global Pd= bus[:, 3]; # Work on it
    global Ybr=zeros(nbr)+zeros(nbr)im;
    global Y=zeros(nbus,nbus)+zeros(nbus,nbus)im;
    aij=ones(nbr);
    for i=1:nbr
        br=network["branch"][string(i)]
        Ybr[i]=1/(br["br_r"]+br["br_x"]im);
        aij[i]=br["tap"]*(cos(br["shift"])+sin(br["shift"]));
        Y[br["f_bus"],br["t_bus"]]=Y[br["f_bus"],br["t_bus"]]-(Ybr[i])/aij[i];
        Y[br["t_bus"],br["f_bus"]]=Y[br["t_bus"],br["f_bus"]]-(Ybr[i])/conj(aij[i]);
        Y[br["f_bus"],br["f_bus"]]=Y[br["f_bus"],br["f_bus"]]+(Ybr[i])/(aij[i]^2);
        Y[br["t_bus"],br["t_bus"]]=Y[br["t_bus"],br["t_bus"]]+(Ybr[i]+br["b_fr"]im);
    end
    global b=imag(Ybr);
    global g=real(Ybr);
    global Ag=zeros(ngen,nbus);
    global Ak=zeros(nbr,nbus);
    global Pmin=zeros(ngen);
    global Pmax=zeros(ngen);
    global Qmin=zeros(ngen);
    global Qmax=zeros(ngen);
    global Vmin=zeros(nbus);
    global Vmax=zeros(nbus);
    global Gen_buses=ones(ngen);
    for i=1:ngen
        Gen_buses[i]=network["gen"][string(i)]["gen_bus"];
        Ag[i,Int32(Gen_buses[i])]=1;
        Pmin[i]=network["gen"][string(i)]["pmin"]*network["gen"][string(i)]["gen_status"]
        Pmax[i]=network["gen"][string(i)]["pmax"]*network["gen"][string(i)]["gen_status"]
        Qmin[i]=network["gen"][string(i)]["qmin"]*network["gen"][string(i)]["gen_status"]
        Qmax[i]=network["gen"][string(i)]["qmax"]*network["gen"][string(i)]["gen_status"]
    end
    global Fmax=zeros(nbr);
    BBr=zeros(nbr,nbr)
    for i=1:nbr
        br=network["branch"][string(i)]
    	Fmax[i]=br["rate_a"]!=0 ? br["rate_a"] : 1000 ;
        Ak[i,br["f_bus"]]=1;
        Ak[i,br["t_bus"]]=-1;
        BBr[i,i]=b[i]
    end
    for i=1:ngen
        gen=network["gen"][string(i)]
        cMin[i]=gen["cost"][1]*Pmin[i]^2+gen["cost"][2]*Pmin[i]+gen["cost"][3];
        if (Pmax[i]>Pmin[i])
            lin_size[i]= (Pmax[i]-Pmin[i])/linearSegments;
            for l=1:linearSegments
                low=Pmin[i]+lin_size[i]*(l-1);
                high=Pmin[i]+lin_size[i]*(l);
                linearCost[i,l]=(gen["cost"][1]*high^2-gen["cost"][1]*low^2)/lin_size[i]+gen["cost"][2];
            end
        end
    end
    global B=imag(Y);
    global G=real(Y);
    global Pd=zeros(nbus);
    global Qd=zeros(nbus);
    global dP=zeros(nbus,nbus);
    global PTDF=zeros(nbr,nbus);
    BT=ones(nbus)
    for j=1:nbus
        bs=network["bus"][string(j)]
        Vmin[j]=bs["vmin"]
        Vmax[j]=bs["vmax"]
        for i=1:length(network["load"])
            Pd[network["load"][string(i)]["load_bus"]]=network["load"][string(i)]["pd"]
            Qd[network["load"][string(i)]["load_bus"]]=network["load"][string(i)]["qd"]
        end
        BT[j]=bs["bus_type"]
    end
    ###_________________________________________
    #TO be added for contingencies
    ###__________________________________________
    #=
    global Slack=findin(BT,3)[1];
    Ak_m=[Ak[:,1:Slack-1] Ak[:,Slack+1:nbus]]
    B_R=[B[1:Slack-1,:]; B[Slack+1:nbus,:]];
    B_C=[B_R[:,1:Slack-1] B_R[:,Slack+1:nbus]]
    B_inv=inv(B_C)
    PTDF_0=zeros(nbr,nbus-1);
    PTDF_0=BBr*Ak_m*B_inv
    PTDF=[PTDF_0[:,1:Slack-1] zeros(nbr,1) PTDF_0[:,Slack:nbus-1]]
    global LODF=zeros(nbr,nbr)
    for k1=1:nbr
        br1=network["branch"][string(k1)]
        for k2=1:nbr
            LODF[k1,k2]=(PTDF[k2,br1["f_bus"]]-PTDF[k2,br1["t_bus"]])/(1-(PTDF[k1,br1["f_bus"]]-PTDF[k1,br1["t_bus"]]))
        end
        LODF[k1,k1]=-1.0
    end=#
end
