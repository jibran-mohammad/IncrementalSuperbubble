#include"superbubbleh1.h"
#include<algorithm>
#include<iostream>

    // read-only view of the current order, minI and maxO (no copy)
    SBView SuperbubbleCalculationH1::view(DFSCONTEXT *sdfs, GraphContext *g) {
        IDFS *d = sdfs->dfsStrategy;
        SBView V;
        V.n = g->displayNumberVertices();
        V.pos = d->reverseTopologicalOrder.data();
        V.at = d->topOrd.data();
        V.minI = d->minI.data();
        V.maxO = d->maxO.data();
        return V;
    }


    void SuperbubbleCalculationH1::computeSB(DFSCONTEXT *sdfs, GraphContext *g) {
        this->tree.build(view(sdfs, g));
        this->pending = false;
    }


    /*
      The graph marks a vertex whose only child is the sink 1 with
      outgoing == 2, and one whose only parent is the source 0 with
      incoming == 2; both flags change when the edge is added, so they are
      read before.
    */
    void SuperbubbleCalculationH1::insertE(DFSCONTEXT *sdfs, GraphContext *g, int u, int v) {
        Graph *G = g->strategy;
        std::map<int,int>::iterator iu = G->reverseMapping.find(u);
        std::map<int,int>::iterator iv = G->reverseMapping.find(v);
        bool uExists = (iu != G->reverseMapping.end());
        bool vExists = (iv != G->reverseMapping.end());
        this->pUWasSink = uExists && G->outgoing[iu->second] == 2;
        this->pVWasSource = vExists && G->incoming[iv->second] == 2;
        this->pDuplicate = false;
        if (uExists && vExists) {
            const std::list<int> &out = G->adjLOutgoing[iu->second];
            this->pDuplicate = (std::find(out.begin(), out.end(), iv->second) != out.end());
        }

    	sdfs->updateEdge(g, u, v);

        this->pU = G->reverseMapping[u];
        this->pV = G->reverseMapping[v];
        this->pCase = sdfs->dfsStrategy->cases;
        this->pending = true;
    }


    void SuperbubbleCalculationH1::updateSB(DFSCONTEXT *sdfs, GraphContext *g) {
        if (!this->pending) return;
        this->pending = false;
        if (this->pDuplicate) {
            this->duplicateEdges++;
            return;
        }
        this->tree.update(view(sdfs, g), this->pCase, this->pU, this->pV, this->pUWasSink, this->pVWasSource);
    }


    std::vector<std::pair<int,int>> SuperbubbleCalculationH1::superbubbles() const {
        std::vector<std::pair<int,int>> out = this->tree.superbubbles();
        out.push_back(std::make_pair(0, 1));
        return out;
    }


    // internal vertex ids are converted to input ids only for printing
    void SuperbubbleCalculationH1::printSuperbubbles(GraphContext *g) {
    	std::cout << "superbubbles of graph are" << std::endl;
    	for (const std::pair<int,int> &sb : this->superbubbles()) {
            std::cout << "<" << g->toOriginal(sb.first) << " " << g->toOriginal(sb.second) << ">" << std::endl;
        }
        std::cout << std::endl;
    }


    // One superbubble per line, "<entrance> <exit>", in input ids.
    void SuperbubbleCalculationH1::writeSuperbubbles(GraphContext *g, std::ostream &out) {
        for (const std::pair<int,int> &sb : this->superbubbles()) {
            out << g->toOriginal(sb.first) << " " << g->toOriginal(sb.second) << "\n";
        }
    }


    int SuperbubbleCalculationH1::numberOfSuperbubbles() {
    	return this->tree.count() + 1;   // + <0 1>
    }
