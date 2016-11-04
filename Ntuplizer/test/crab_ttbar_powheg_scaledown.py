from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.transferOutputs = True
config.General.requestName = 'TTBar_Powheg_ScaleDown_v7'

config.section_('JobType')
config.JobType.psetName = 'hadronic.py'
config.JobType.pluginName = 'Analysis'
config.JobType.inputFiles = ['JEC/START53_V27_L1FastJet_AK7PFchs.txt', 'JEC/START53_V27_L2Relative_AK7PFchs.txt', 'JEC/START53_V27_L3Absolute_AK7PFchs.txt', 'JEC/START53_V27_Uncertainty_AK7PFchs.txt',
                                'JEC/START53_V27_L1FastJet_AK5PFchs.txt', 'JEC/START53_V27_L2Relative_AK5PFchs.txt', 'JEC/START53_V27_L3Absolute_AK5PFchs.txt', 'JEC/START53_V27_Uncertainty_AK5PFchs.txt',
                                'JEC/Winter14_V5_DATA_L1FastJet_AK7PFchs.txt', 'JEC/Winter14_V5_DATA_L2Relative_AK7PFchs.txt', 'JEC/Winter14_V5_DATA_L3Absolute_AK7PFchs.txt',
                                'JEC/Winter14_V5_DATA_L2L3Residual_AK7PFchs.txt', 'JEC/Winter14_V5_DATA_Uncertainty_AK7PFchs.txt',
                                'JEC/Winter14_V5_DATA_L1FastJet_AK5PFchs.txt', 'JEC/Winter14_V5_DATA_L2Relative_AK5PFchs.txt', 'JEC/Winter14_V5_DATA_L3Absolute_AK5PFchs.txt',
                                'JEC/Winter14_V5_DATA_L2L3Residual_AK5PFchs.txt', 'JEC/Winter14_V5_DATA_Uncertainty_AK5PFchs.txt']
config.JobType.pyCfgParams = ['runOnData=0', 'runOnCrab=1']

config.section_('Data')
config.Data.inputDataset = '/TT_scaledown_CT10_TuneZ2star_8TeV-powheg-tauola/asady-TTJets_scaledown-0454285497e05eceaa5af276c1e8e4ef/USER'
config.Data.unitsPerJob = 5
config.Data.splitting = 'FileBased'
config.Data.outputDatasetTag = 'TTBar_Powheg_ScaleDown_v7'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T3_US_FNALLPC'