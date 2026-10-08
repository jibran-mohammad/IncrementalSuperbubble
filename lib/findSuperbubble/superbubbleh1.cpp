#include"superbubbleh1.h"

    void SuperbubbleCalculationH1::computeSB(DFSCONTEXT *sdfs,GraphContext *g){
	    int t1,t;
	    std::stack<ppii> S;
	    for(t1=0;t1<=g->displayNumberVertices()-1;t1++){
	    	t=sdfs->dfsStrategy->topOrd[t1];
	    	/*
	    	if(sdfs->dfsStrategy->minI[t] == INT_MAX) {
	    		while(!S.empty()){
	    			S.pop();
	    		}
	    	}
	    	*/
	    	// if minI[t] goes beyond what is on top of the stack, then pop elements until the elemement is found
	    	// which does where minI[t] is within the bounds.
		    while(!S.empty() && sdfs->dfsStrategy->reverseTopologicalOrder[S.top().first] > sdfs->dfsStrategy->reverseTopologicalOrder[sdfs->dfsStrategy->minI[t]]){
			    ppii s1 = S.top();
			    S.pop();
			    this->nodeList[S.top().first].push_back(s1.first);
	            if(!S.empty()) {
			        ppii s2 = S.top();
			        S.pop();
			        if(s1.second.first == INT_MAX) {
			        	S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(s1.second.first, sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
			        else {
			        S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.first], sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
			    }

		    }
			
			//Superbubble condition
		    if(!S.empty() && S.top().second.second == t){
			    //assert(MIN(S.top().second.first,sdfs->dfsStrategy->minI[t])==S.top().first);
			    //this->SB.push_back(pii(S.top().first,t));
			    //printf("<%d,%d>\n", S.top().first,t);
			    SBID id_st = nextSBID++;
			    SB[id_st] = pii(S.top().first,t);
                sbAtStart[S.top().first] = id_st;

                SuperbubbleNode node(id_st);
                int maxChildLevel = 0;
                
                // adding node into the tree with parent children relationship
                for (auto [si, ti] : innerSB[S.top().first]) {
                    SBID childID = sbAtStart[si];

                    SBTree[childID].parent = id_st;
                    SBTree[id_st].children.push_back(childID);

                    maxChildLevel = std::max(maxChildLevel, SBTree[childID].level);
                }

                SBTree[id_st].level = maxChildLevel + 1;
                sbOfVertex[S.top().first] = id_st;
                sbOfVertex[t] = id_st;  
                if (!innerSB[s].empty()) {
                    // has inner superbubbles → assign only free vertices
                    for (int v : nodeList[s]) {
                        if (sbOfVertex[v] == INVALID_SBD) {
                            sbOfVertex[v] = id_st;
                        }
                    }
                } else {
                    // leaf superbubble → assign all internal vertices
                    for (int v : nodeList[s]) {
                        sbOfVertex[v] = id_st;
                    }
                }
			    ppii s1 = S.top();
			    S.pop();
			    this->innerSB[S.top().first].push_back(pii(s1.top().first, t))
			    this->nodeList[S.top().first].push_back(s1.first);
			    if(!S.empty()){
				    ppii s2 = S.top();
				    S.pop();
			        if(s1.second.first == INT_MAX) {
			        	S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(s1.second.first, sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
			        else {
			            S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.first], sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
				       
			       
			    }

			    
		    }
		    

		    // minI value change
		    if(!S.empty()){
			    ppii s1 = S.top();
			    S.pop();
		        if(s1.second.first == INT_MAX) {
		        	s1.second.first = sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[sdfs->dfsStrategy->minI[t]],s1.second.first)];
		        }
		        else {
			        s1.second.first = sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[sdfs->dfsStrategy->minI[t]],sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.first])];
			    }
			    S.push(s1);
		    }
		
		    S.push(ppii(t,pii(INT_MAX,sdfs->dfsStrategy->maxO[t])));
	    }	

    }

    
    void SuperbubbleCalculationH1::printSuperbubbles() {
    	std::cout << "superbubbles of graph are" << std::endl;
    	for (pii element : this->SB) {
            std::cout << "<"<<element.first <<" " << element.second <<">";
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }


    int SuperbubbleCalculationH1::numberOfSuperbubbles() {
    	int numElements = this->SB.size();
    	return numElements;
    }

    
    void SuperbubbleCalculationH1::insertE(DFSCONTEXT *sdfs,GraphContext *g,int u, int v) {
    	SBID rootId;
    	sdfs->updateEdge(g,u,v);
        rootId = this->findLCA(sbOfVertex[u], sbOfVertex[v]);
        this->clearSubtree(rootId);
        this->recomputeInside(sdfs,g,u,v);

    }

    
    SBID SuperbubbleCalculationH1::findLCA(SBID u, SBID v) {
        // equalize levels
        while (SBTree[u].level > SBTree[v].level)
            u = SBTree[u].parent;

        while (SBTree[v].level > SBTree[u].level)
            v = SBTree[v].parent;

        // climb together
        while (u != v) {
            u = SBTree[u].parent;
            v = SBTree[v].parent;
        }

        return u;
    }

    void clearSubtree(SBID rootID) {
        std::vector<SBID> subtree;
        this->collectSubtree(rootID, subtree);

        for (SBID id : subtree) {
            int s = SB[id].first;    

            // clear sbOfVertex links
            for (int v : nodeList[s]) {
                sbOfVertex[v] = INVALID_SBD;
            }

            // clear construction buffers
            nodeList[s].clear();
            innerSB[s].clear();

            // clear start mapping
            sbAtStart[s] = INVALID_SBD;
        }

        // erase SBs and tree nodes
        for (SBID id : subtree) {
            SBTree.erase(id);
            SB.erase(id);
        }
    
    }

    void collectSubtree(SBID rootID, std::vector<SBID> &subtree) {
    	std::stack<SBID> st;
        st.push(rootID);

        while (!st.empty()) {
            SBID x = st.top(); st.pop();
            subtree.push_back(x);

            for (SBID c : SBTree[x].children) {
                st.push(c);
            }
        }
    }

    void SuperbubbleCalculationH1::recomputeInside(DFSCONTEXT *sdfs,GraphContext *g, int u, int v) {
         int t1,t;
	    std::stack<ppii> S;
	    for(t1=sdfs->dfsStrategy->reverseTopologicalOrder[u];t1<=sdfs->dfsStrategy->reverseTopologicalOrder[v];t1++){
	    	t=sdfs->dfsStrategy->topOrd[t1];
	    	/*
	    	if(sdfs->dfsStrategy->minI[t] == INT_MAX) {
	    		while(!S.empty()){
	    			S.pop();
	    		}
	    	}
	    	*/
	    	// if minI[t] goes beyond what is on top of the stack, then pop elements until the elemement is found
	    	// which does where minI[t] is within the bounds.
		    while(!S.empty() && sdfs->dfsStrategy->reverseTopologicalOrder[S.top().first] > sdfs->dfsStrategy->reverseTopologicalOrder[sdfs->dfsStrategy->minI[t]]){
			    ppii s1 = S.top();
			    S.pop();
			    this->nodeList[S.top().first].push_back(s1.first);
	            if(!S.empty()) {
			        ppii s2 = S.top();
			        S.pop();
			        if(s1.second.first == INT_MAX) {
			        	S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(s1.second.first, sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
			        else {
			        S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.first], sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
			    }

		    }
			
			//Superbubble condition
		    if(!S.empty() && S.top().second.second == t){
			    //assert(MIN(S.top().second.first,sdfs->dfsStrategy->minI[t])==S.top().first);
			    this->SB.push_back(pii(S.top().first,t));
			    //printf("<%d,%d>\n", S.top().first,t);
			    SBID id_st = nextSBID++;
			    SB[id_st] = pii(S.top().first,t);
                sbAtStart[S.top().first] = id_st;

                SuperbubbleNode node(id_st);
                int maxChildLevel = 0;
                
                // adding node into the tree with parent children relationship
                for (auto [si, ti] : innerSB[S.top().first]) {
                    SBID childID = sbAtStart[si];

                    SBTree[childID].parent = id_st;
                    SBTree[id_st].children.push_back(childID);

                    maxChildLevel = std::max(maxChildLevel, SBTree[childID].level);
                }

                SBTree[id_st].level = maxChildLevel + 1;
                sbOfVertex[S.top().first] = id_st;
                sbOfVertex[t] = id_st  
                if (!innerSB[s].empty()) {
                    // has inner superbubbles → assign only free vertices
                    for (int v : nodeList[s]) {
                        if (sbOfVertex[v] == INVALID_SBD) {
                            sbOfVertex[v] = id_st;
                        }
                    }
                } else {
                    // leaf superbubble → assign all internal vertices
                    for (int v : nodeList[s]) {
                        sbOfVertex[v] = id_st;
                    }
                }
			    ppii s1 = S.top();
			    S.pop();
			    this->innerSB[S.top().first].push_back(pii(s1.top().first, t))
			    this->nodeList[S.top().first].push_back(s1.first);
			    if(!S.empty()){
				    ppii s2 = S.top();
				    S.pop();
			        if(s1.second.first == INT_MAX) {
			        	S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(s1.second.first, sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
			        else {
			            S.push(ppii(s2.first,pii(sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.first], sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.first])],
			            sdfs->dfsStrategy->topOrd[MAX(sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.second],sdfs->dfsStrategy->reverseTopologicalOrder[s2.second.second])]  )));
			        }
				       
			       
			    }

			    
		    }
		    

		    // minI value change
		    if(!S.empty()){
			    ppii s1 = S.top();
			    S.pop();
		        if(s1.second.first == INT_MAX) {
		        	s1.second.first = sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[sdfs->dfsStrategy->minI[t]],s1.second.first)];
		        }
		        else {
			        s1.second.first = sdfs->dfsStrategy->topOrd[MIN(sdfs->dfsStrategy->reverseTopologicalOrder[sdfs->dfsStrategy->minI[t]],sdfs->dfsStrategy->reverseTopologicalOrder[s1.second.first])];
			    }
			    S.push(s1);
		    }
		
		    S.push(ppii(t,pii(INT_MAX,sdfs->dfsStrategy->maxO[t])));
	    }	
    }




