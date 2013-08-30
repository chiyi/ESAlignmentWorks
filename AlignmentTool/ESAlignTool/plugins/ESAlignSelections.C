bool ESAlignTool::pass_TrackSelection(reco::TrackCollection::const_iterator itTrack)
//bool ESAlignTool::pass_TrackSelection(reco::Track *itTrack)
{
 bool res=false;
 if(   itTrack->pt()>1.5
     &&fabs(itTrack->outerZ())>260&&fabs(itTrack->outerZ())<280
     &&fabs(itTrack->eta())<2.3&&fabs(itTrack->eta())>1.7
     &&itTrack->found()>=10 //NHits
     && ((itTrack->qualityMask())%8) >= 4 //highPurity Bit
   )
  res=true;
 return res;
}   

