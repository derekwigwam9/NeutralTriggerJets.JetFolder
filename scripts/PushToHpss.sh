# 'PushToHpss.sh'
#
# Use this to copy data to HPSS.

date="19Oct2016"

echo "Copying to HPSS! Today is $date"

htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/Pythia17.unfolding."$date".tar ./output/Pythia17*.root >& Pythia17.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/Pythia19.unfolding."$date".tar ./output/Pythia19*.root >& Pythia19.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/PriorTestBay.ChargePt."$date".tar ./output/PriorTestBay/ChargePt* >& PriorTestBay.ChargePt.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/PriorTestBay.test."$date".tar ./output/PriorTestBay/test.* >& PriorTestBay.test.hpss.log &
htar -r /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/PriorTestBay.test."$date".tar ./output/PriorTestBay/ExpBay.Sep19.log >& PriorTestBay.test.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/PriorTestSvd.ChargePt."$date".tar ./output/PriorTestSvd/ChargePt* >& PriorTestSvd.ChargePt.hpss.log &
htar -r /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/PriorTestSvd.ChargePt."$date".tar ./output/PriorTestSvd/ChargeTest.* >& PriorTestSvd.ChargePt.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/PriorTestSvd.test."$date".tar ./output/PriorTestSvd/test.* >& PriorTestSvd.test.hpss.log &
htar -r /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/PriorTestSvd.test."$date".tar ./output/PriorTestSvd/ExpSVD.Sep21.log >& PriorTestSvd.test.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/pp200aug.Pythia17g.unfolding."$date".tar ./output/Pythia17g_pp200aug/pp200aug* >& pp200aug.Pythia17g.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/pp200aug.Pythia17p.unfolding."$date".tar ./output/Pythia17p_pp200aug/pp200aug* >& pp200aug.Pythia17p.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/pp200aug.Pythia19g.unfolding."$date".tar ./output/Pythia19g_pp200aug/pp200aug* >& pp200aug.Pythia19g.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/pp200aug.Pythia19p.unfolding."$date".tar ./output/Pythia19p_pp200aug/pp200aug* >& pp200aug.Pythia19p.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/pp200aug.Pythia20g.unfolding."$date".tar ./output/Pythia20g_pp200aug/pp200aug* >& pp200aug.Pythia20g.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/ResponseMatrices/Pythia17.response."$date".tar ./output/ResponseMatrices/Pythia17r.* >& Pythia17r.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/ResponseMatrices/Pythia19.response."$date".tar ./output/ResponseMatrices/Pythia19r.* >& Pythia19r.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/ResponseMatrices/Pythia20.response."$date".tar ./output/ResponseMatrices/Pythia20r.* >& Pythia20r.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/ResponseMatrices/Pythia22.response."$date".tar ./output/ResponseMatrices/Pythia22r.* >& Pythia22r.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/ResponseMatrices/debug."$date".tar ./output/ResponseMatrices/Debug.root >& debug.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/UnfoldingTests."$date".tar ./output/tests/*.root >& UnfoldingTests.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/ToyUnfolding."$date".tar ./output/ToyUnfolding/*.root >& ToyUnfolding.hpss.log &
htar -cf /home/d/dmawxc/NeutralTriggeredJets/JetUnfolding/input."$date".tar ./input/*.root >& input.hpss.log &
