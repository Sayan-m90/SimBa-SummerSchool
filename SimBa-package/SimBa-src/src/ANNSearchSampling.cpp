 #include "ANNSearchSampling.h" 

//namespace DASComp = AbstractSimplicialComplex;
namespace ANNSearch
{
	// randomly select a point from a point set
	int RandomSelect(const std::set<int>& intSet)
	{
		
		int r;
		int size = (int)intSet.size();
		r = rand();
		r = r % size;
		std::set<int>::const_iterator csIter = intSet.begin();
		for (int i = 0; i < r; i++)
			csIter++;
		return *csIter;
	} 
	//
	void Subsampling_EuclideanDistance( const double sqDelta,
										const PointSet &ptSet,
										std::map<int, int> &SubsampleIndices,
										std::vector<int> &ColorMapping
										)
	{//
		//
		int nPts = (int)ptSet._PointSet->size();
		int dim = ptSet._dimension;
		//
		int k = 0; // parameter k in k-ann 
		//
		double eps = 0.0; // distance tolerance
		//
		ANNdist sqRadius = sqDelta;
		//
		SubsampleIndices.clear(); // initial empty subsampling
		//
		std::vector<double> curDistToSubsampling(nPts);
		//
		//initialize color mapping with negative values to indicate no mapping yet
		//
		memset(&(ColorMapping[0]), -1, sizeof(int) * nPts);
		//
		std::unordered_set<int> uncolored_points;
		//
		std::vector<int> vecRandomSequence;
		for (int i = 0; i < nPts; i++)
		{
			uncolored_points.insert(i);
			vecRandomSequence.push_back(i);
		}
		//generating random sequence
		for (int i = 0; i < nPts; ++i)
		{
			int size = nPts - i;
			int j = rand();
			j = j % size + i;
			//swap
			int temp = vecRandomSequence[i];
			vecRandomSequence[i] = vecRandomSequence[j];
			vecRandomSequence[j] = temp;
		}
		int curPos = 0;
		// 
		int subsampling_index = 0;
		/*************/
		ANNpoint queryPt;
		ANNidxArray nnIdx = NULL;
		ANNdistArray dists = NULL;

		queryPt = annAllocPt(dim); // need to delloc
		
 		nnIdx = new ANNidx[nPts];
		dists = new ANNdist[nPts];
		/**************************/

		ANNpointArray dataPts;
		dataPts = annAllocPts(nPts, dim);
		// adding vertices
 
		for (int i = 0; i < nPts; i++){
			for (int j = 0; j < dim; j++)
			{
				dataPts[i][j] = (*ptSet._PointSet)[i][j]; 
			}
		}
		///////////////////
		ANNkd_tree* kdTree;
		kdTree = new ANNkd_tree(dataPts, 
								nPts,
								dim);
		//
		int NeighborCounter = 0;
		//
		// perform the sub-sampling here
		while (!uncolored_points.empty())
		{
			int cur_pt_idx /*= RandomSelect(uncolored_points)*/;
			for (; curPos < nPts; ++curPos)
			{
				if (uncolored_points.find(vecRandomSequence[curPos]) != uncolored_points.end())
				{
					cur_pt_idx = vecRandomSequence[curPos];
					curPos += 1;
					break;
				}
			}
			//
			// set up the query point as this point
			for (int j = 0; j < dim; j++)
			{
				queryPt[j] = (*ptSet._PointSet)[cur_pt_idx][j];
			}
			// get the number of points inside the fixed radius
			k = kdTree->annkFRSearch(
									queryPt,
									sqRadius,
									0, // to get the number of points within the radius
									nnIdx,
									dists,
									eps);
			NeighborCounter = kdTree->annkFRSearch(
											queryPt,
											sqRadius,
											k,
											nnIdx,
											dists,
											eps);
			// record edges
			if (NeighborCounter > k)
			{
				std::cout << "k is too small" << std::endl;
				exit(0);
			}
			if (NeighborCounter > 0)
			{// add edges
				if (nnIdx[0] != cur_pt_idx && dists[0] > ZERO)
				{
					std::cout << "SELF is NOT contained" << std::endl;
					std::cout << "i " << cur_pt_idx << std::endl;
					std::cout << "nn " << nnIdx[0] << std::endl;
					std::cout << dists[0] << std::endl;
					exit(0);
				}
				else
				{
					ColorMapping[cur_pt_idx] = subsampling_index;
					uncolored_points.erase(cur_pt_idx);
				}
				for (int j = 1; j < NeighborCounter; j++)
				{
					if (ColorMapping[nnIdx[j]] < 0)
					{ // this point is not colored
						// color this point
						ColorMapping[nnIdx[j]] = subsampling_index;
						// remove it from uncolored vertex set
						uncolored_points.erase(nnIdx[j]);
						// record its distance to sub-sampling
						curDistToSubsampling[nnIdx[j]] = dists[j];
					}
					else
					{// this point is colored before
					// update its color based the distance
						if (dists[j] < curDistToSubsampling[nnIdx[j]] )
							// by default , if dists[j] == curDistToSubsampling[nnIdx[j]], we assign it with smaller index point
						{
							curDistToSubsampling[nnIdx[j]] = dists[j];
							ColorMapping[nnIdx[j]] = subsampling_index;
						}
					}
				}
			}
			else
			{
				std::cout << "isolated point !!!" << std::endl;
				exit(0);
			}
			//
			SubsampleIndices[cur_pt_idx] = subsampling_index++;
			//
			//std::cout << "i " << i << std::endl;
		}
 		// clean data
		//annDeallocPts(dataPts);

		delete [] nnIdx;
		delete [] dists;
		delete kdTree;
		if (dataPts)
			annDeallocPts(dataPts);
		if (queryPt)
			annDeallocPt(queryPt);
		annClose();
	}
											
	//
	void BuildDistanceGraph_ann( SimpleGraph& retGraph, 
								const PointSet &ptSet,
								const double sqDistance)
	{
		int nPts = (int)ptSet._PointSet->size();
		int dim =  ptSet._dimension;
		//
		int k = 2000;//nPts;
		if (k > nPts)
			k = nPts;
		double eps = 0.0;
		ANNdist sqRad = sqDistance;
		//
		retGraph.InitNodes(nPts);

		/*************/
		ANNpoint queryPt;
		ANNidxArray nnIdx = NULL;
		ANNdistArray dists = NULL;

		queryPt = annAllocPt(dim); // need to delloc
		
 		nnIdx = new ANNidx[nPts];
		dists = new ANNdist[nPts];
		/**************************/

		ANNpointArray dataPts;
		dataPts = annAllocPts(nPts, dim);
		// adding vertices
 
		for (int i = 0; i < nPts; i++){
			for (int j = 0; j < dim; j++)
			{
				dataPts[i][j] = (*ptSet._PointSet)[i][j]; 
			}
		}
		///////////////////
		ANNkd_tree* kdTree;
		kdTree = new ANNkd_tree(dataPts, 
								nPts,
								dim);
		//
		int NeighborCounter = 0;
		for (int i = 0; i < nPts; i++)
		{
			for (int j = 0; j < dim; j++)
			{
				queryPt[j] = (*ptSet._PointSet)[i][j];
			}
		//
			// get the number of points inside the fixed radius
			k = kdTree->annkFRSearch(
									queryPt,
									sqRad,
									0, // to get the number of points within the radius
									nnIdx,
									dists);//,
									//eps); omitting eps para means exact search
			NeighborCounter = kdTree->annkFRSearch(
											queryPt,
											sqRad,
											k,
											nnIdx,
											dists,
											eps);
			// record edges
			if (NeighborCounter > k)
			{
				std::cout << "k is too small" << std::endl;
				exit(0);
			}
			if (NeighborCounter > 0)
			{// add edges
				if (nnIdx[0] != i)
				{
					std::cout << "SELF is NOT contained" << std::endl;
					std::cout << "i " << i << std::endl;
					std::cout << "nn " << nnIdx[0] << std::endl;
					std::cout << dists[0] << std::endl;
					exit(0);
				}
				for (int j = 1; j < NeighborCounter; j++)
				{
					//if (!retGraph.vecNode[i].is_neighbor(nnIdx[j]))
					if (!retGraph.is_edge(i, nnIdx[j]))
					{
						retGraph.AddEdge(i, nnIdx[j]);
 					}
				}
			}
			//std::cout << "i " << i << std::endl;
		}
 		// clean data
		//annDeallocPts(dataPts);

		delete [] nnIdx;
		delete [] dists;
		delete kdTree;
		if (dataPts)
			annDeallocPts(dataPts);
		if (queryPt)
			annDeallocPt(queryPt);
		annClose();
	} 

/************************************************/
 
//
	void SetColorMappingAndExtractColoredGraph(	const std::vector<int>& ColorMapping,
												const SimpleGraph& RipsGraph,
 												SimpleGraph& colorGraph							
							)
 	{
 		// iterate over the edges in RipsGraph
		//
		colorGraph.InitNodes((int)RipsGraph.vecNode.size());
		//
		for (int i = 0; i < (int) RipsGraph.vecNode.size(); i++)
		{
			// set color mapping
			colorGraph.vecNode[i].color = ColorMapping[i];
			//RipsGraph.vecNode[i].color = ColorMapping[i];
			//
			if (!RipsGraph.vecNode[i].edgeWeights.empty())
			{
				for (std::map<int, float>::const_iterator mIter = RipsGraph.vecNode[i].edgeWeights.begin();
					mIter != RipsGraph.vecNode[i].edgeWeights.end();
					mIter++)
				{
					if (ColorMapping[mIter->first] != ColorMapping[i])
					{ // bi-colored edge
						colorGraph.AddEdge(i, mIter->first);
					}
				}
			}
		} // for i size
	}

};
