// IMPORTS SECTION
import React, { useState, useEffect, Suspense, lazy } from 'react';
import { 
  Table, 
  Card,
  Form,
  Input,
  InputNumber,
  Space,
  Modal,
  Typography,
  Layout,
  Button,
  message,
  Spin,
  Select,
  TreeSelect,
  Row,
  Col,
  Tooltip,
  Slider,
  Statistic,
  Checkbox,
  Divider,
  Image,
  Tag,
  Menu,
  Anchor,
  Tabs
} from 'antd';
import { 
  SearchOutlined, 
  DownloadOutlined, 
  ClearOutlined, 
  EyeOutlined, 
  DeploymentUnitOutlined,
  BulbOutlined,
  FilterOutlined,
  FileTextOutlined,
  HomeOutlined,
  QuestionCircleOutlined,
  InfoCircleOutlined,
  BookOutlined,
  DesktopOutlined
} from '@ant-design/icons';
import { Routes, Route, Link, useLocation } from 'react-router-dom';
import axios from 'axios';
import { 
  BarChart, 
  Bar, 
  XAxis, 
  YAxis, 
  CartesianGrid, 
  Tooltip as ChartTooltip, 
  ResponsiveContainer, 
  Cell,
  PieChart,
  Pie,
  Legend,
  ScatterChart,
  Scatter,
  ZAxis
} from 'recharts';

const MolstarViewer = lazy(() => import('./MolstarViewer'));

const { Header, Content, Footer } = Layout;
const { Title, Text, Paragraph } = Typography;
const { Option } = Select;

import { API_BASE } from './types';

const CHART_COLORS = [
  '#13c2c2', '#52c41a', '#1890ff', '#722ed1', '#eb2f96',
  '#faad14', '#f5222d', '#fa8c16', '#a0d911', '#1d39c4',
  '#08979c', '#d4380d', '#531dab', '#c41d7f', '#d48806'
];

// MENU SECTION
const HelpPage: React.FC = () => (
  <Row justify="center">
    <Col xs={24} lg={20} xl={16}>
      <Card style={{ borderRadius: '8px', boxShadow: '0 2px 8px rgba(0,0,0,0.05)' }}>

        <div id="intro" style={{ marginBottom: 60 }}>
          <Title level={2} style={{ textAlign: 'center' }}>RNAbridge User Guide</Title>
          <Paragraph style={{ fontSize: '16px', textAlign: 'center', color: '#666' }}>
            RNAbridge is an advanced platform for identifying and analyzing RNA extended, non-canonical helices. 
            This guide will help you fully utilize the tool's capabilities.
          </Paragraph>
        </div>

        <Divider />

        {/* SECTION 1: Search */}
        <div id="search" style={{ marginBottom: 60 }}>
          <Title level={4}><SearchOutlined /> 1. Search & Filtering Panel</Title>
          <Paragraph>
            The search panel provides multiple ways to narrow down the dataset:
          </Paragraph>
          <ul style={{ marginBottom: 24 }}>
            <li><strong>CIF ID:</strong> Search for specific structures using their 4-character identifier.</li>
            <li><strong>Sequence Pattern:</strong> Find motifs containing specific nucleotide sequences (supports A, C, G, U).</li>
            <li><strong>Type:</strong> Filter by complexity, ranging from simple helices to multi-way junctions (up to 8+ way).</li>
            <li><strong>Coaxially Stacked:</strong> Filter junctions by identifying specific pairs of stems (e.g., Stem 1 and Stem 3) that exhibit coaxial stacking.</li>
            <li><strong>Angle (0-50°):</strong> Filter structures based on their global curvature, which often indicates functional strain.</li>
            <li><strong>NT Range:</strong> Limit results by the total number of nucleotides within the identified motif.</li>
          </ul>
          <Image 
            src="/guide/search_panel.png" 
            alt="Search panel"
            style={{ borderRadius: '8px', border: '1px solid #f0f0f0', marginBottom: 16 }}
            fallback="https://placehold.co/1000x300?text=Search+Panel+Screenshot"
          />
          <Paragraph type="secondary">
            Tip: Filtering is fully automatic for most inputs. The results table updates instantly as you adjust any parameter, except for the sequence pattern which requires pressing enter or clicking the search icon.
          </Paragraph>
        </div>

        <Divider />

        {/* SECTION 2: Statistics */}
        <div id="stats" style={{ marginBottom: 60 }}>
          <Title level={4}><BulbOutlined /> 2. Statistics Dashboard</Title>
          <Paragraph>
            The dashboard provides a high-level overview of your current search results, allowing you to quickly assess the distribution and composition of the identified motifs:
          </Paragraph>

          <Row gutter={[16, 16]}>
            <Col xs={24} md={8}>
              <Card title="Results Count" size="small" type="inner" style={{ height: '100%' }}>
                <Paragraph>Shows the total number of records matching your active filters.</Paragraph>
                <Image 
                  src="/guide/stats_count.png" 
                  alt="Results count"
                  style={{ borderRadius: '4px' }}
                  fallback="https://placehold.co/300x150?text=Count+Stat"
                />
              </Card>
            </Col>
            <Col xs={24} md={16}>
              <Card title="Bend Angle Distribution" size="small" type="inner" style={{ height: '100%' }}>
                <Paragraph>A bar chart showing how motifs are distributed across different bend angle ranges (0° to 50°). <strong>Hover over any bar</strong> to see the exact count, or <strong>click a bar</strong> to instantly filter the results to that specific angle range.</Paragraph>
                <Image 
                  src="/guide/stats_histogram.png" 
                  alt="Angle distribution"
                  style={{ borderRadius: '4px' }}
                  fallback="https://placehold.co/600x150?text=Histogram+Chart"
                />
              </Card>
            </Col>
            <Col span={24}>
              <Card title="Composition Analysis" size="small" type="inner">
                <Paragraph>A pie chart representing the proportion of different motif types (Helices vs. various Junctions) in your selection. <strong>Hover over a slice</strong> to view the precise count for that specific motif type.</Paragraph>
                <div style={{ textAlign: 'center' }}>
                  <Image 
                    src="/guide/stats_pie.png" 
                    alt="Composition pie chart"
                    style={{ borderRadius: '4px', maxHeight: '200px' }}
                    fallback="https://placehold.co/800x200?text=Pie+Chart+Screenshot"
                  />
                </div>
              </Card>
            </Col>
          </Row>
          <Paragraph type="secondary" style={{ marginTop: 16 }}>
            Note: These charts are interactive and update automatically whenever filters are applied.
          </Paragraph>
        </div>

        <Divider />

        {/* SECTION 3: Records */}
        <div id="table" style={{ marginBottom: 60 }}>
          <Title level={4}><FileTextOutlined /> 3. Result Analysis</Title>
          <Paragraph>
            The results table displays all identified motifs matching your search criteria. Each row provides a summary of the structural data:
          </Paragraph>
          <ul style={{ marginBottom: 24 }}>
            <li><strong>CIF ID:</strong> The 4-character PDB entry identifier (e.g., <Tag color="blue">8XTP</Tag>). Clicking it opens the original RCSB PDB page.</li>
            <li><strong>Source / Molecule:</strong> Displays the name of the molecule and the organism from which the RNA was derived.</li>
            <li><strong>Method:</strong> The experimental technique used to determine the structure (e.g., X-ray Diffraction, Cryo-EM).</li>
            <li><strong>Res. (Å):</strong> The experimental resolution of the structure in Angstroms. Lower values indicate higher detail.</li>
            <li><strong>Nts Count:</strong> Total number of nucleotides contained within the identified motif.</li>
            <li><strong>Bend Angle (°):</strong> The global curvature of the motif. High angles indicate sharp structural bends.</li>
            <li><strong>Type:</strong> Categorizes the motif as either a <strong>HELIX</strong> or a <strong>JUNCTION</strong>, along with its specific segment count.</li>
            <li><strong>Visualization:</strong> Quick access to the <strong>Preview</strong> button, which opens the integrated 2D/3D structure viewer.</li>
            <li><strong>Download:</strong> Each row contains direct links to download the <strong>CIF</strong> (3D coordinates) and <strong>PML</strong> (PyMOL visualization script) files.</li>
          </ul>

          <div style={{ marginBottom: 32 }}>
            <Text strong>Record Example (3-way Junction):</Text>
            <Image 
              src="/guide/table_row.png" 
              alt="Single record example"
              style={{ borderRadius: '8px', border: '1px solid #f0f0f0', marginTop: 8 }}
              fallback="https://placehold.co/1000x150?text=8XTP+3-way-junction+Record+Screenshot"
            />
          </div>

          <Divider dashed />

          {/* SUB-SECTION: Expansion Details */}
          <div id="details" style={{ marginTop: 32 }}>
            <Title level={5}>Detailed Motif Analysis (Expansion View)</Title>
            <Paragraph>
              Clicking anywhere on a row expands the record to reveal in-depth structural parameters. You can also download the raw data by clicking the <strong>Download JSON</strong> button in the top-right corner of the expansion.
            </Paragraph>

            <Row gutter={[16, 16]}>
              <Col xs={24} md={12}>
                <Card title="Helix Components" size="small" type="inner" style={{ height: '100%' }}>
                  <Paragraph>Breakdown of individual segments within the helix:</Paragraph>
                  <ul style={{ fontSize: '13px', paddingLeft: '20px' }}>
                    <li><strong>Type:</strong> Identifies the part (Stem, Hairpin, etc.). Bulges are further classified as <Tag color="orange">BULGE-IN</Tag> or <Tag color="orange">BULGE-OUT</Tag>.</li>
                    <li><strong>Size (2D / 3D):</strong> Compares schematic length (2D) vs. actual nucleotides in the stacking path (3D).</li>
                    <li><strong>Local Angle:</strong> The curvature specific to this segment.</li>
                    <li><strong>Stacking Path:</strong> Displays the specific nucleotide sequence involved in coaxial stacking (e.g., 101→102→105).</li>
                    <li><strong>Sequence:</strong> The primary nucleotide sequence of the segment.</li>
                  </ul>
                  <Image 
                    src="/guide/details_helix.png" 
                    alt="Helix details"
                    style={{ borderRadius: '4px', marginTop: 8 }}
                    fallback="https://placehold.co/500x250?text=Helix+Details+Screenshot"
                  />
                </Card>
              </Col>
              <Col xs={24} md={12}>
                <Card title="Coaxial Stacking (Junctions)" size="small" type="inner" style={{ height: '100%' }}>
                  <Paragraph>Analysis of stem interactions within the junction:</Paragraph>
                  <ul style={{ fontSize: '13px', paddingLeft: '20px' }}>
                    <li><strong>Pair ID:</strong> Identifies the two stems involved in the coaxial interaction (e.g., STEM 1 - STEM 3).</li>
                    <li><strong>Chain:</strong> The molecular chain identifier.</li>
                    <li><strong>Angle:</strong> The relative geometric orientation angle between the two stacked stems.</li>
                    <li><strong>Sequence 1 / 2:</strong> A detailed 2D visual representation of the stacked stems, including nucleotide numbering and base-pairing markers.</li>
                  </ul>
                  <Image 
                    src="/guide/details_junction.png" 
                    alt="Junction details"
                    style={{ borderRadius: '4px', marginTop: 8 }}
                    fallback="https://placehold.co/500x250?text=Junction+Details+Screenshot"
                  />
                </Card>
              </Col>
            </Row>
          </div>
        </div>

        <Divider />

        {/* SECTION 4: Visualization */}
        <div id="viewer" style={{ marginBottom: 60 }}>
          <Title level={4}><EyeOutlined /> 4. Integrated Structural Preview (3D & 2D)</Title>
          <Paragraph>
            The <strong>Structure Preview</strong> tool provides a synchronized side-by-side view of the identified motif, bridging secondary topology with tertiary geometry.
          </Paragraph>

          <Card size="small" type="inner" style={{ borderRadius: '8px' }}>
            <ul style={{ marginBottom: 16 }}>
              <li><strong>2D Schematic:</strong> An interactive diagram displaying base-pairing and connectivity. Loop core strands and motifs are color-coded using distinct vibrant palettes.</li>
              <li><strong>3D Model (Mol*):</strong> A full spatial representation, perfectly synchronized with the 2D view. Components use the same colors, allowing for instant identification of structural features.</li>
              <li><strong>Component Legend:</strong> A comprehensive legend displayed below the viewers, showing both the Helix and Junction palettes simultaneously for easy reference.</li>
            </ul>

            <div style={{ textAlign: 'center', background: '#f5f5f5', borderRadius: '4px', padding: '10px' }}>
              <Image 
                src="/guide/preview_combined.png" 
                alt="Combined 2D/3D Preview"
                style={{ width: '100%', borderRadius: '4px', border: '1px solid #ddd' }}
                fallback="https://placehold.co/1000x600?text=Combined+2D+and+3D+Preview+Screenshot"
              />
            </div>
          </Card>
        </div>

        <Divider />

        {/* SECTION 5: Export */}
        <div id="export" style={{ marginBottom: 40 }}>
          <Title level={4}><DownloadOutlined /> 5. Data Export & Batch Downloads</Title>
          <Paragraph>
            RNAbridge supports multiple ways to save your findings for offline research:
          </Paragraph>
          <ul style={{ marginBottom: 24 }}>
            <li><strong>Export to CSV:</strong> Generates a spreadsheet containing the metadata and geometric parameters for all currently filtered records.</li>
            <li><strong>Download All Results (.zip):</strong> A powerful feature that packages the <strong>CIF coordinates</strong> and <strong>PML scripts</strong> for <em>every</em> result in your current search. This allows you to instantly open and analyze multiple structures in PyMOL.</li>
          </ul>
          <Image 
            src="/guide/export_options.png" 
            alt="Export options"
            style={{ borderRadius: '8px', border: '1px solid #f0f0f0', marginBottom: 16 }}
            fallback="https://placehold.co/1000x100?text=Export+Options+Screenshot"
          />
        </div>

      </Card>
    </Col>
  </Row>
);

const AboutPage: React.FC = () => (
  <Row justify="center">
    <Col xs={24} lg={20} xl={16}>
      <Card style={{ borderRadius: '8px', boxShadow: '0 2px 8px rgba(0,0,0,0.05)' }}>
        <Title level={2} style={{ textAlign: 'center', marginBottom: '40px' }}>About RNAbridge</Title>
        
        <section style={{ marginBottom: '40px' }}>
          <Title level={4}>Mission & Motivation</Title>
          <Paragraph style={{ fontSize: '16px', lineHeight: '1.7' }}>
            RNAbridge is a specialized repository and analysis platform designed to bridge the gap between complex 3D RNA structures and actionable geometric insights. While databases like the Protein Data Bank (PDB) provide raw coordinates, researchers often require specific parameters—such as <strong>global and local bend angles</strong>, <strong>segment orientations</strong>, and <strong>coaxial stacking patterns</strong>—to fully understand RNA folding and biological function.
          </Paragraph>
          <Paragraph style={{ fontSize: '16px', lineHeight: '1.7' }}>
            Our mission is to provide an accessible "bridge" to this structural data, offering a pre-analyzed collection of RNA helices and junctions extracted from the entire PDB and processed with high-precision bioinformatics tools.
          </Paragraph>
        </section>

        <Divider />

        <section style={{ marginBottom: '40px' }}>
          <Title level={4}>Technical Pipeline</Title>
          <Paragraph>
            The RNAbridge database is built on a robust automated pipeline:
          </Paragraph>
          <ul style={{ lineHeight: '2' }}>
            <li><strong>Data Acquisition:</strong> Automated extraction and filtering of RNA-containing structures from the Protein Data Bank.</li>
            <li><strong>Motif Identification:</strong> Systematic detection of RNA helices and multi-way junctions (up to 8-way and beyond).</li>
            <li><strong>Geometric Characterization:</strong> Calculation of local and global bend angles, integrating specialized analysis scripts and the <em>x3dna</em> toolkit.</li>
            <li><strong>Stacking Analysis:</strong> Automated detection of coaxial stacking patterns to identify structural continuity between stems in complex junctions.</li>
          </ul>
        </section>

        <Divider />

        <section style={{ marginBottom: '40px' }}>
          <Title level={4}>Research Team</Title>
          <ul style={{ lineHeight: '2' }}>
            <li><a href="https://www.linkedin.com/in/damian-zakrzewski2004/" target="_blank" rel="noopener noreferrer" style={{ color: 'inherit' }}>Damian Zakrzewski</a></li>
            <li>Tomasz Zok</li>
            <li>Maciej Antczak</li>
          </ul>
        </section>

        <Divider />

        <section style={{ textAlign: 'center', padding: '20px 0' }}>
          <Title level={4}>Affiliations</Title>
          <div style={{ textAlign: 'center' }}>
            <Text strong style={{ fontSize: '16px' }}>Poznan University of Technology</Text><br />
            <Text type="secondary">Institute of Computing Science</Text>
          </div>
        </section>

      </Card>
    </Col>
  </Row>
);

const CitePage: React.FC = () => (
  <Card title={<span><BookOutlined /> Cite Us</span>} style={{ borderRadius: '8px' }}>
    <Typography>
      <Paragraph>Citation information will be available soon.</Paragraph>
    </Typography>
  </Card>
);


interface Result {
  id: number;
  pdb_id: string;
  molecule?: string;
  organism?: string;
  resolution?: number;
  method?: string;
  total_nt: number;
  global_bend_angle: number | null;
  segment_count_folder: string;
  path_svg: string | null;
  path_cif: string;
  path_pml: string;
  path_json: string;
  type: 'helix' | 'junction';
  details?: any;
}

const App: React.FC = () => {
  const location = useLocation();
  const [form] = Form.useForm();
  const [results, setResults] = useState<Result[]>([]);
  const [totalResults, setTotalResults] = useState<number>(0);
  const [searchStats, setSearchStats] = useState<any>(null);
  const [selectedRowKeys, setSelectedRowKeys] = useState<React.Key[]>([]);
  const [availableIds, setAvailableIds] = useState<string[]>([]);
  const [maxNtLimit, setMaxNtLimit] = useState<number>(500);
  const [loading, setLoading] = useState<boolean>(false);
  const [currentPage, setCurrentPage] = useState<number>(1);
  const [pageSize, setPageSize] = useState<number>(10);
  const [sortState, setSortState] = useState<{ field: string | null, order: 'asc' | 'desc' | null }>({ field: null, order: null });
  const [availableTypes, setAvailableTypes] = useState<{value: string, label: string}[]>([]);
  const [currentSvg, setCurrentSvg] = useState<string | null>(null);
  const [svgPreviewVisible, setSvgPreviewVisible] = useState<boolean>(false);

  const [maxWayFromStats, setMaxWayFromStats] = useState<number>(0);
  const selectedTypes = Form.useWatch('segment_type', form);

  const stemOptions = React.useMemo(() => {
    let maxToDisplay = maxWayFromStats;

    if (selectedTypes && selectedTypes.length > 0) {
      const counts = selectedTypes.map((t: string) => {
        const match = t.match(/^(\d+)/);
        return match ? parseInt(match[1], 10) : 0;
      });
      maxToDisplay = Math.max(...counts);
    }

    return Array.from({ length: maxToDisplay }, (_, i) => ({
      value: `STEM_${i + 1}`,
      label: `STEM ${i + 1}`
    }));
  }, [selectedTypes, maxWayFromStats]);

  // Clear stem selection if it becomes invalid after changing type
  useEffect(() => {
    const s1 = form.getFieldValue('stacking_stem1');
    const s2 = form.getFieldValue('stacking_stem2');
    let changed = false;
    const newVals: any = {};

    if (s1 && !stemOptions.find(o => o.value === s1)) {
      newVals.stacking_stem1 = undefined;
      changed = true;
    }
    if (s2 && !stemOptions.find(o => o.value === s2)) {
      newVals.stacking_stem2 = undefined;
      changed = true;
    }

    if (changed) {
      form.setFieldsValue(newVals);
    }
  }, [stemOptions, form]);

  const [combinedModalVisible, setCombinedModalVisible] = useState<boolean>(false);
  const [currentRecord, setCurrentRecord] = useState<Result | null>(null);
  const [highlightData, setHighlightData] = useState<any>(null);

  const openCombinedModal = (record: Result) => {
    setCurrentRecord(record);
    
    // Logic for highlights
    const components = record.details?.components;
    if (components && Array.isArray(components)) {
      const stackingIds: string[] = [];
      if (record.type === 'junction' && Array.isArray(record.details?.coaxial_pairs)) {
        record.details.coaxial_pairs.flat().forEach((s: any) => {
          stackingIds.push(s.toString().toLowerCase());
        });
      }

      const highlightsList = components.flatMap((c: any) => {
        const residues = c.residues || [];
        const cid = c.id?.toLowerCase() || '';
        const isStacking = stackingIds.includes(cid);
        
        let finalColor = { r: 180, g: 180, b: 180 };
        let label = c.id || 'Segment';
        let priority = 0;

        if (c.color) {
          finalColor = c.color;
          priority = 2; 
        } else if (isStacking) {
          finalColor = { r: 0, g: 0, b: 0 };
          label = `Coaxial Stem: ${c.id}`;
          priority = 1;
        }

        return residues.map((r: any) => ({
          start: r.start,
          end: r.end,
          chain: r.chain || c.chain || 'A',
          start_icode: r.start_icode,
          end_icode: r.end_icode,
          color: finalColor,
          label: label,
          priority: priority
        }));
      });
      setHighlightData(highlightsList.sort((a: any, b: any) => a.priority - b.priority));
    } else {
      setHighlightData(null);
    }
    
    setCombinedModalVisible(true);
  };

  const [jsonPreviewVisible, setJsonPreviewVisible] = useState<boolean>(false);
  const [currentJsonContent, setCurrentJsonContent] = useState<string>('');

  const [csvPreviewVisible, setCsvPreviewVisible] = useState<boolean>(false);
  const [csvContent, setCsvContent] = useState<string>('');

  const onSearch = async (values: any, page = 1, limit = 10, sorter: any = null) => {
    setLoading(true);
    setCurrentPage(page);
    setPageSize(limit);

    let currentSortField = sortState.field;
    let currentSortOrder = sortState.order;

    if (sorter && sorter.order) {
      currentSortField = sorter.field;
      currentSortOrder = sorter.order === 'descend' ? 'desc' : 'asc';
      setSortState({ field: currentSortField, order: currentSortOrder });
    } else if (sorter && !sorter.order && sorter !== 'reset') {
      // Obsługa wyczyszczenia sortowania (trzecie kliknięcie)
      currentSortField = null;
      currentSortOrder = null;
      setSortState({ field: null, order: null });
    } else if (sorter === 'reset') {
      currentSortField = null;
      currentSortOrder = null;
      setSortState({ field: null, order: null });
    }
	    
    try {
      const { angle_range, nt_range, ...rest } = values;

      // Use URLSearchParams to handle multiple values for the same key (e.g. ?segment_type=a&segment_type=b)
      const params = new URLSearchParams();

      if (rest.search_pdb) params.append('search_pdb', rest.search_pdb);
      if (rest.sequence) params.append('sequence', rest.sequence);
      if (rest.stacking_stem1) params.append('stacking_stem1', rest.stacking_stem1);
      if (rest.stacking_stem2) params.append('stacking_stem2', rest.stacking_stem2);

      if (rest.segment_type && Array.isArray(rest.segment_type)) {
        rest.segment_type.forEach((t: string) => params.append('segment_type', t));
      }

      params.append('min_angle', (angle_range ? angle_range[0] : 0).toString());
      params.append('max_angle', (angle_range ? angle_range[1] : 50).toString());
      params.append('min_nt', (nt_range ? nt_range[0] : 0).toString());
      params.append('max_nt', (nt_range ? nt_range[1] : maxNtLimit).toString());

      params.append('page', page.toString());
      params.append('limit', limit.toString());

      if (currentSortField) {
        params.append('sort_by', currentSortField);
        params.append('sort_order', currentSortOrder || 'asc');
      }

      const response = await axios.get(`${API_BASE}/api/search`, { params });      
      if (response.data && response.data.results) {
        setResults(response.data.results);
        setTotalResults(response.data.total);
        setSearchStats(response.data.stats); // Zapisujemy statystyki z serwera
        if (response.data.results.length === 0 && page === 1) {
          message.warning('No results found for the selected criteria.');
        }
      } else {
        const data = Array.isArray(response.data) ? response.data : [];
        setResults(data);
        setTotalResults(data.length);
      }
    } catch (error) {
      console.error('Error fetching data:', error);
      message.error('Failed to connect to the API server.');
    } finally {
      setLoading(false);
    }
  };

  const handleBatchDownload = async () => {
    setLoading(true);
    try {
      const values = form.getFieldsValue();
      const { angle_range, nt_range, ...rest } = values;
      const params = new URLSearchParams();
      if (rest.search_pdb) params.append('search_pdb', rest.search_pdb);
      if (rest.sequence) params.append('sequence', rest.sequence);
      if (rest.stacking_stem1) params.append('stacking_stem1', rest.stacking_stem1);
      if (rest.stacking_stem2) params.append('stacking_stem2', rest.stacking_stem2);
      if (rest.segment_type && Array.isArray(rest.segment_type)) {
        rest.segment_type.forEach((t: string) => params.append('segment_type', t));
      }
      params.append('min_angle', (angle_range ? angle_range[0] : 0).toString());
      params.append('max_angle', (angle_range ? angle_range[1] : 50).toString());
      params.append('min_nt', (nt_range ? nt_range[0] : 0).toString());
      params.append('max_nt', (nt_range ? nt_range[1] : maxNtLimit).toString());

      const response = await axios.get(`${API_BASE}/api/export-zip`, { 
        params,
        responseType: 'blob'
      });

      const url = window.URL.createObjectURL(new Blob([response.data]));
      const link = document.createElement('a');
      link.href = url;
      link.setAttribute('download', `RNAbridge_all_results_${new Date().toISOString().slice(0, 10)}.zip`);
      document.body.appendChild(link);
      link.click();
      link.remove();
      message.success(`Successfully downloaded all filtered structures.`);
    } catch (error) {
      console.error('Download error:', error);
      message.error('Failed to generate ZIP archive.');
    } finally {
      setLoading(false);
    }
  };

  const handleBarClick = (data: any) => {
    if (!data || !data.range) return; 
    const matches = data.range.match(/\d+/g);
    if (!matches || matches.length < 2) return;
    const min = Number(matches[0]);
    let max = Number(matches[1]);
    if (data.range.endsWith(')')) {
      max -= 0.01;
    }
    const currentValues = form.getFieldsValue();
    const newFilters = { ...currentValues, angle_range: [min, max] };
    form.setFieldsValue({ angle_range: [min, max] });
    onSearch(newFilters);
    message.info(`Filtered by angle range: ${data.range}`);
  };

  const handlePieClick = (data: any) => {
    if (!data || !data.name) return;
    const groupName = data.name;
    const typesToFilter = getTypesForGroup(groupName); 
    const currentValues = form.getFieldsValue();
    const newFilters = { ...currentValues, segment_type: typesToFilter };
    form.setFieldsValue({ segment_type: typesToFilter });
    onSearch(newFilters);
    message.info(`Filtered by type: ${groupName}`);
  };

  const exportToCsv = async () => {
    setLoading(true);
    try {
      const values = form.getFieldsValue();
      const { angle_range, nt_range, ...rest } = values;
      const params = new URLSearchParams();
      if (rest.search_pdb) params.append('search_pdb', rest.search_pdb);
      if (rest.sequence) params.append('sequence', rest.sequence);
      if (rest.stacking_stem1) params.append('stacking_stem1', rest.stacking_stem1);
      if (rest.stacking_stem2) params.append('stacking_stem2', rest.stacking_stem2);
      if (rest.segment_type && Array.isArray(rest.segment_type)) {
        rest.segment_type.forEach((t: string) => params.append('segment_type', t));
      }
      params.append('min_angle', (angle_range ? angle_range[0] : 0).toString());
      params.append('max_angle', (angle_range ? angle_range[1] : 50).toString());
      params.append('min_nt', (nt_range ? nt_range[0] : 0).toString());
      params.append('max_nt', (nt_range ? nt_range[1] : maxNtLimit).toString());

      if (sortState.field) {
        params.append('sort_by', sortState.field);
        params.append('sort_order', sortState.order || 'asc');
      }

      const response = await axios.get(`${API_BASE}/api/export-csv`, { params });
      setCsvContent(response.data);
      setCsvPreviewVisible(true);
    } catch (error) {
      console.error('CSV Export error:', error);
      message.error('Failed to export CSV data.');
    } finally {
      setLoading(false);
    }
  };

  const downloadCsv = () => {
    const blob = new Blob([csvContent], { type: 'text/csv;charset=utf-8;' });
    const link = document.createElement('a');
    const url = URL.createObjectURL(blob);
    link.setAttribute('href', url);
    link.setAttribute('download', `RNAbridge_export_${new Date().toISOString().slice(0, 10)}.csv`);
    link.style.visibility = 'hidden';
    document.body.appendChild(link);
    link.click();
    document.body.removeChild(link);
    setCsvPreviewVisible(false);
    message.success('CSV file downloaded successfully.');
  };

  const fetchIds = async () => {
    try {
      const response = await axios.get(`${API_BASE}/api/ids`);
      setAvailableIds(response.data);
    } catch (error) {
      console.error('Error fetching IDs:', error);
    }
  };

  const fetchStats = async () => {
    try {
      const response = await axios.get(`${API_BASE}/api/stats`);
      const maxVal = response.data.max_nt;
      const maxWay = response.data.max_way || 0;
      const folders = response.data.folders || []; 

      setMaxNtLimit(maxVal);
      form.setFieldsValue({ nt_range: [0, maxVal] });
      setMaxWayFromStats(maxWay);
    
      const formatLabel = (folder: string) => {
        if (folder === '8plus-junctions') return '8+ way Junctions';
        
        const segMatch = folder.match(/^(\d+)-segment/);
        if (segMatch) {
          const num = parseInt(segMatch[1], 10);
          return `Helices (${num} segment${num > 1 ? 's' : ''})`;
        }

        const wayMatch = folder.match(/^(\d+)-way/);
        if (wayMatch) {
          return `${wayMatch[1]}-way Junctions`;
        }

        return folder; 
      };

      const types = folders
        .map((f: string) => ({
          value: f,
          label: formatLabel(f),
          isHelix: f.includes('segment'),
          num: parseInt(f.match(/\d+/)?.[0] || '0', 10)
        }))
        .sort((a, b) => {
          if (a.isHelix && !b.isHelix) return -1;
          if (!a.isHelix && b.isHelix) return 1;
          return a.num - b.num;
        })
        .map(item => ({ value: item.value, label: item.label }));
        
      setAvailableTypes(types);  
    } catch (error) {
      console.error('Error fetching stats:', error);
    }
  }; 
  const getGroupForType = (folderName: string) => {
    if (!folderName || folderName.includes('1-segment')) return null;

    if (folderName.includes('segment')) {
      const match = folderName.match(/^(\d+)/);
      if (match) {
        const num = parseInt(match[1], 10);
        if (num === 2) return '2 seg';
        if (num === 3) return '3 seg';
        if (num >= 4) return '4 seg+';
      }
    }

    if (folderName.includes('way') || folderName.includes('plus')) {
      const match = folderName.match(/^(\d+)/);
      if (match) {
        const num = parseInt(match[1], 10);
        if (num === 3) return '3 way';
        if (num === 4) return '4 way';
        if (num === 5) return '5 way';
        if (num >= 6) return '6 way+';
      }
      if (folderName.includes('8plus') || folderName.includes('way-junction')) {
        return '6 way+';
      }
    }
    return null;
  };

  const getTypesForGroup = (groupName: string) => {
    return availableTypes
      .map(t => t.value)
      .filter(val => getGroupForType(val) === groupName);
  };

  const prepareHistogramData = () => {
    const typesGroups = ['2 seg', '3 seg', '4 seg+', '3 way', '4 way', '5 way', '6 way+'];
    const bins = Array.from({ length: 10 }, (_, i) => {
      const min = i * 5;
      const max = (i + 1) * 5;
      const rangeLabel = i === 9 ? `[${min}, ${max}]` : `[${min}, ${max})`;
      
      return {
        range: rangeLabel,
        ...Object.fromEntries(typesGroups.map(t => [t, 0]))
      };
    });

    if (!searchStats || !searchStats.angles) return bins;

    searchStats.angles.forEach((item: any) => {
      const binIndex = Math.min(Math.floor(item.bin / 5), 9);
      // Używamy naszego tłumacza zamiast sztywnej listy!
      const groupName = getGroupForType(item.folder || "");
      
      if (groupName && binIndex >= 0 && binIndex < 10) {
        bins[binIndex][groupName] = (bins[binIndex][groupName] as number) + item.count;
      }
    });
    return bins;
  };

  const preparePieData = () => {
    const counts: Record<string, number> = {
      '2 seg': 0, '3 seg': 0, '4 seg+': 0,
      '3 way': 0, '4 way': 0, '5 way': 0, '6 way+': 0
    };

    if (!searchStats || !searchStats.pie) return [];

    Object.entries(searchStats.pie).forEach(([f, count]: [string, any]) => {
      const groupName = getGroupForType(f);
      if (groupName) {
        counts[groupName] += count;
      }
    });

    return Object.entries(counts)
      .filter(([_, val]) => val > 0)
      .map(([name, value]) => ({ name, value }));
  };
  
  const TYPE_COLORS: Record<string, string> = {
    '2 seg': '#13c2c2',
    '3 seg': '#52c41a',
    '4 seg+': '#237804',
    '3 way': '#1890ff',
    '4 way': '#722ed1',
    '5 way': '#eb2f96',
    '6 way+': '#faad14'
  };

  const PIE_COLORS = ['#13c2c2', '#52c41a', '#237804', '#1890ff', '#722ed1', '#eb2f96', '#faad14'];

  useEffect(() => {
    fetchIds();
    fetchStats();
    onSearch({});
  }, []);

  const clearFilters = () => {
    form.resetFields();
    onSearch({}, 1, 10, 'reset');
  };

  const runExample = (params: any) => {
    form.setFieldsValue(params);
    onSearch(params);
  };

  const openSvgPreview = (svgPath: string) => {
    setCurrentSvg(`${API_BASE}/api/files/${svgPath}`);
    setSvgPreviewVisible(true);
  };

  const openJsonPreview = async (jsonPath: string) => {
    try {
      const response = await axios.get(`${API_BASE}/api/files/${jsonPath}`);
      setCurrentJsonContent(JSON.stringify(response.data, null, 2));
      setJsonPreviewVisible(true);
    } catch (error) {
      console.error('Error loading JSON:', error);
      message.error('Failed to load JSON file.');
    }
  };

  const columns = [
    {
      title: 'CIF ID',
      dataIndex: 'pdb_id',
      key: 'pdb_id',
      render: (text: string) => (
        <a href={`https://www.rcsb.org/structure/${text}`} target="_blank" rel="noopener noreferrer">
          <Tag color="blue" style={{ fontWeight: 'bold', cursor: 'pointer' }}>{text}</Tag>
        </a>
      ),
      sorter: true,
    },
    {
      title: 'Source / Molecule',
      key: 'molecule_info',
      width: 250,
      render: (_: any, record: Result) => (
        <Space direction="vertical" size={0}>
          <a href={`https://www.rcsb.org/structure/${record.pdb_id}`} target="_blank" rel="noopener noreferrer">
            <Text strong style={{ fontSize: '12px', color: '#1890ff', cursor: 'pointer' }}>
              {record.molecule || 'N/A'}
            </Text>
          </a>
          <Text type="secondary" style={{ fontSize: '11px', fontStyle: 'italic' }}>{record.organism || 'Unknown Organism'}</Text>
        </Space>
      ),
    },
    {
      title: 'Method',
      dataIndex: 'method',
      key: 'method',
      width: 120,
      render: (text: string) => <Text style={{ fontSize: '11px' }}>{text || 'N/A'}</Text>,
    },
    {
      title: 'Res. (Å)',
      dataIndex: 'resolution',
      key: 'resolution',
      width: 80,
      render: (val: number | null) => val ? <Text code>{val.toFixed(2)}</Text> : <Text type="secondary">-</Text>,
      sorter: true,
    },
    {
      title: 'Nts Count',
      dataIndex: 'total_nt',
      key: 'total_nt',
      render: (val: number) => <Text strong>{val}</Text>,
      sorter: true,
    },
    {
      title: 'Bend Angle (°)',
      dataIndex: 'global_bend_angle',
      key: 'global_bend_angle',
      render: (val: number | null, record: Result) => {
        if (val === 0 && record.segment_count_folder?.includes('1-segment')) {
          return '-';
        } 
        return val !== null ? <Tag color={val > 45 ? 'volcano' : 'cyan'}>{val.toFixed(2)}°</Tag> : '-';
      },
      sorter: true,
    },
    {
      title: 'Type',
      dataIndex: 'segment_count_folder',
      key: 'segment_count_folder',
      render: (text: string, record: Result) => (
        <Space direction="vertical" size={0}>
          <Tag color={record.type === 'helix' ? 'green' : 'orange'}>{record.type.toUpperCase()}</Tag>
          <Text type="secondary" style={{ fontSize: '11px' }}>{text}</Text>
        </Space>
      ),
    },

    {
      title: 'Visualization',
      key: 'viz',
      width: 120,
      render: (_: any, record: Result) => (
        <Button
          type="primary"
          onClick={(e) => { e.stopPropagation(); openCombinedModal(record); }}
          style={{ width: '100%', height: '36px', fontWeight: 'bold' }}
        >
          Preview
        </Button>
      ),
    },
    {
      title: 'Download',
      key: 'files',
      render: (_: any, record: Result) => (
        <Space direction="vertical" onClick={(e) => e.stopPropagation()}>
          <Button 
            size="small" 
            icon={<DownloadOutlined />} 
            href={`${API_BASE}/api/files/${record.path_cif}`} 
            download
            style={{ width: '80px' }}
          >
            CIF
          </Button>
          <Button 
            size="small" 
            icon={<DownloadOutlined />} 
            href={`${API_BASE}/api/files/${record.path_pml}`} 
            download
            style={{ width: '80px' }}
          >
            PML
          </Button>
        </Space>
      ),
    },
  ];

  const expandedRowRender = (record: Result) => {
    if (record.type === 'helix') {
      const segments = record.details?.db_segments || [];
      const segmentColumns = [
        { title: 'ID', dataIndex: 's_id', key: 's_id', width: 60 },
        {
          title: 'Type',
          dataIndex: 'type',
          key: 'type',
          render: (t: string, seg: any) => {
            const isStem = t.toLowerCase().includes('stem');
            if (t.toUpperCase() === 'BULGE') {
              const s = seg.stacking;
              if (s === 'BULGE-IN' || s === 'BULGE-OUT') {
                return <Tag color="orange">{s}</Tag>;
              }
            }
            return <Tag color={isStem ? 'blue' : 'orange'}>{t}</Tag>;
          }
        },
        { 
          title: 'Size (2D / 3D)', 
          key: 'size_display',
          render: (_: any, seg: any) => {
            const type = seg.type.toLowerCase();
            if (type.includes('stem')) {
              return <Text type="secondary" style={{ fontSize: '12px' }}>{seg.size_2d}</Text>;
            }
            if (type.includes('hairpin')) {
              const size = seg.size_2d.split('x')[0];
              return <Tag color="blue" style={{ margin: 0 }}>{size} nt</Tag>;
            }
            return (
              <Space size={4}>
                <Tag color="blue" style={{ margin: 0, minWidth: '45px', textAlign: 'center' }}>
                  <b>{seg.size_2d}</b>
                </Tag>
                <Tag color="green" style={{ margin: 0, minWidth: '45px', textAlign: 'center' }}>
                  <b>{seg.size_3d}</b>
                </Tag>
              </Space>
            );
          }
        },
        { 
          title: 'Local Angle', 
          dataIndex: 'bend_angle', 
          key: 'bend_angle', 
          width: 100,
          render: (v: number, record: any) => {
            const isStem = record.type.toLowerCase().includes('stem');
            if (isStem && (v === 0 || v === null)) return '-';
            return v !== null ? v.toFixed(1) + '°' : '-';
          }
        },
        { 
          title: 'Stacking Path', 
          dataIndex: 'stacking_path', 
          key: 'stacking_path',
          width: 140,
          render: (path: any, seg: any) => {
            if (!path || !Array.isArray(path) || path.length === 0) return '-';
            
            const formatPath = (p: any[]) => p.join('→');
            
            let displayStr = '';
            let fullTooltip = '';

            if (Array.isArray(path[0])) {
              // Internal loop (multiple strands)
              return (
                <Tooltip title={path.map((p: any[], i: number) => `${i + 1}: ${formatPath(p)}`).join('  ')}>
                  <div style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                    {path.map((p: any[], i: number) => (
                      <Text key={i} code style={{ fontSize: '11px', cursor: 'help', whiteSpace: 'nowrap' }}>
                        {p.length > 0 ? `${i + 1}: ${p[0]}...${p[p.length-1]}` : `${i+1}: []`}
                      </Text>
                    ))}
                  </div>
                </Tooltip>
              );
            } else {
              // Bulge / Hairpin
              displayStr = path.length > 3 ? `${path[0]}→${path[1]}...${path[path.length-1]}` : formatPath(path);
              fullTooltip = formatPath(path);
              return (
                <Tooltip title={fullTooltip}>
                  <Text code style={{ fontSize: '11px', cursor: 'help' }}>
                    {displayStr}
                  </Text>
                </Tooltip>
              );
            }
          }
        },
        { 
          title: 'Sequence', 
          key: 'sequence_display',
          render: (_: any, seg: any) => (
            <Text style={{ fontFamily: 'monospace', fontWeight: 'bold', fontSize: '15px', letterSpacing: '1px', color: '#001529' }}>
              {seg.sequence}
            </Text>
          )
        },
      ];

      return (
        <div style={{ margin: '10px 0', padding: '15px', background: '#fafafa', borderRadius: '4px', border: '1px solid #eee' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '10px' }}>
            <Title level={5} style={{ margin: 0 }}>Helix Components</Title>            <Button 
              size="small" 
              icon={<DownloadOutlined />} 
              href={`${API_BASE}/api/files/${record.path_json}`} 
              download
            >
              Download JSON
            </Button>
          </div>
          <Table 
            columns={segmentColumns} 
            dataSource={segments} 
            pagination={false} 
            size="small" 
            rowKey={(s, idx) => `segment-${idx}`}
            bordered
          />
        </div>
      );
    } else {
      const details = record.details || {};
      const pairs = details.coaxial_pairs || [];
      const angles = details.all_angles || {};
      const stackingPaths = details.stacking_paths || {};
      const stems = details.stems || {};


      const dataSource = (Array.isArray(pairs) ? pairs : []).map((pair: any, idx: number) => {
        if (!pair || pair.length < 2) return null;
        const s1_key = pair[0].toLowerCase();
        const s2_key = pair[1].toLowerCase();
        
        const stem1 = stems[s1_key];
        const stem2 = stems[s2_key];
        
        const id1 = s1_key.split('_').pop();
        const id2 = s2_key.split('_').pop();
        const angle = angles[`stem_${id1}_${id2}`] ?? angles[`stem_${id2}_${id1}`];
	const path = stackingPaths[`stem_${id1}_${id2}`] ?? stackingPaths[`stem_${id2}_${id1}`];
        return {
          key: `stack-pair-${idx}`,
          id: `${pair[0].toUpperCase()} - ${pair[1].toUpperCase()}`,
          chain: stem1?.strand5p?.first?.chain || '-',
          angle: angle,
	  path,
          stem1,
          stem2
        };
      }).filter(Boolean);

      const renderSingleBridge = (stem: any) => {
        if (!stem) return null;
        const s5 = stem.strand5p || {};
        const s3 = stem.strand3p || {};
        const seq5 = s5.sequence || '';
        const seq3 = s3.sequence || '';
        const seq3Rev = seq3.split('').reverse().join('');

        const n5_f = s5.first?.number ?? '?';
        const n5_l = s5.last?.number ?? '?';
        const n3_f = s3.first?.number ?? '?';
        const n3_l = s3.last?.number ?? '?';

        const prefix5 = `${n5_f}-`;
        const prefix3 = `${n3_l}-`;
        
        // Obliczamy przesunięcie, aby kreski były nad literami
        const p5_len = prefix5.length;
        const p3_len = prefix3.length;
        const max_p = Math.max(p5_len, p3_len);

        const pad5 = ' '.repeat(max_p - p5_len);
        const pad3 = ' '.repeat(max_p - p3_len);
        const pipePad = ' '.repeat(max_p);
        
        const pipes = '|'.repeat(seq5.length);

        return (
          <div style={{ 
            fontFamily: 'monospace', 
            fontSize: '14px', 
            fontWeight: 'bold', 
            lineHeight: '1.1',
            whiteSpace: 'pre',
            color: '#001529'
          }}>
            {pad5}{prefix5}{seq5}-{n5_l}<br />
            {pipePad}{pipes}<br />
            {pad3}{prefix3}{seq3Rev}-{n3_f}
          </div>
        );
      };

      const junctionColumns = [
        { 
          title: 'Pair ID', 
          dataIndex: 'id', 
          key: 'id', 
          width: 180,
          render: (text: string) => <Text strong>{text}</Text>
        },
        { title: 'Chain', dataIndex: 'chain', key: 'chain', width: 80 },
        { 
          title: 'Angle', 
          dataIndex: 'angle', 
          key: 'angle',
          width: 100,
          render: (v: any) => (typeof v === 'number') ? (
            <Tag color={v <= 50 ? 'green' : 'cyan'}>{v.toFixed(1)}°</Tag>
          ) : '-'
        },
        { 
          title: 'Sequence 1', 
          key: 'bridge1',
          render: (_: any, record: any) => renderSingleBridge(record.stem1)
        },
	{ 
          title: 'Sequence 2', 
          key: 'bridge2',
          render: (_: any, record: any) => renderSingleBridge(record.stem2)
        },
        { 
          title: 'Stacking Path', 
          dataIndex: 'path', 
          key: 'path',
          width: 140,
          render: (path: any) => {
            if (!path || !Array.isArray(path) || path.length === 0) return '-';

            const formatPath = (p: any[]) => p.join('→');

            if (Array.isArray(path[0])) {
              // Multiple independent stacking paths between the same stem pair
              return (
                <Tooltip title={path.map((p: any[], i: number) => `${i + 1}: ${formatPath(p)}`).join('  ')}>
                  <div style={{ display: 'flex', flexDirection: 'column', gap: '2px' }}>
                    {path.map((p: any[], i: number) => (
                      <Text key={i} code style={{ fontSize: '11px', cursor: 'help', whiteSpace: 'nowrap' }}>
                        {p.length > 0 ? `${i + 1}: ${p[0]}...${p[p.length-1]}` : `${i+1}: []`}
                      </Text>
                    ))}
                  </div>
                </Tooltip>
              );
            } else {
              const displayStr = path.length > 3 ? `${path[0]}→${path[1]}...${path[path.length-1]}` : formatPath(path);
              return (
                <Tooltip title={formatPath(path)}>
                  <Text code style={{ fontSize: '11px', cursor: 'help' }}>
                    {displayStr}
                  </Text>
                </Tooltip>
              );
            }
          }
        },
      ];
      return (
        <div style={{ margin: '10px 0', padding: '15px', background: '#fafafa', borderRadius: '4px', border: '1px solid #eee' }}>
          <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '15px' }}>
            <Title level={5} style={{ margin: 0 }}>Coaxial Stacking Analysis</Title>
            <Button size="small" icon={<DownloadOutlined />} href={`${API_BASE}/api/files/${record.path_json}`} download>
              Download JSON
            </Button>
          </div>
          <Table
            columns={junctionColumns}
            dataSource={dataSource}
            pagination={false}
            size="small"
            rowKey="key"
            bordered
            locale={{ emptyText: 'No coaxial stacking detected' }}
          />
        </div>
      );
    }
  };

  return (
    <Layout style={{ minHeight: '100vh', backgroundColor: '#f0f2f5' }}>
      <Header style={{ background: '#001529', padding: '0 20px', display: 'flex', alignItems: 'center' }}>
        <Space style={{ marginRight: '40px' }}>
          <img src="/rnabridge.svg" alt="RNAbridge Logo" style={{ height: '50px', verticalAlign: 'middle' }} />
        </Space>
        <Menu 
          theme="dark" 
          mode="horizontal" 
          selectedKeys={[location.pathname]} 
          style={{ flex: 1, minWidth: 0, borderBottom: 'none' }}
          items={[
            { key: '/', icon: <HomeOutlined />, label: <Link to="/">Home</Link> },
            { key: '/help', icon: <QuestionCircleOutlined />, label: <Link to="/help">Help</Link> },
            { key: '/about', icon: <InfoCircleOutlined />, label: <Link to="/about">About</Link> },
            { key: '/cite', icon: <BookOutlined />, label: <Link to="/cite">Cite Us</Link> },
          ]}
        />
      </Header>

      <Content style={{ padding: '15px 15px' }}>
        <Routes>
          <Route path="/" element={
            <Row gutter={[16, 16]}>
              <Col span={24}>
                <Card style={{ borderRadius: '8px' }}>
                  <Title level={4}>About RNAbridge</Title>
                  <Text style={{ fontSize: '16px' }}>
                    RNAbridge is a specialized platform designed for the identification and analysis of RNA structural motifs, focusing on helices and junctions. It enables researchers to explore complex RNA architectures, identify coaxial stacking patterns, and access detailed geometric data for various RNA segments to support structural biology research.
                  </Text>
                </Card>
              </Col>

              <Col span={24}>
                <Card
                  title={<span><BulbOutlined style={{ color: '#faad14' }} /> Quick Examples</span>}
                  style={{ borderRadius: '8px' }}
                  styles={{ body: { padding: '15px 20px' } }}
                >
                  <Space wrap>
                    <Button type="dashed" onClick={() => runExample({ angle_range: [40, 50] })}>
                      Example 1: Sharp Bends (40-50°)
                    </Button>
                    <Button type="dashed" onClick={() => runExample({ segment_type: ['4-way-junctions'], nt_range: [70, maxNtLimit], angle_range: [0, 20] })}>
                      Example 2: Large 4-way Junctions (NT &gt; 70, 0-20°)
                    </Button>
                    <Button type="dashed" onClick={() => runExample({ segment_type: ['8plus-junctions'] })}>
                      Example 3: Complex Junctions (8+ segments)
                    </Button>
                    <Button type="dashed" onClick={() => runExample({ segment_type: ['4-segment-helis'], nt_range: [0, 50], angle_range: [25, 50] })}>
                      Example 4: Helices (4 seg), Angle 25-50°, NT &le; 50
                    </Button>
                  </Space>
                </Card>
              </Col>

          <Col span={24}>
            <Card
              title={<span><FilterOutlined /> Filters</span>}
              style={{ borderRadius: '8px' }}
              styles={{ body: { padding: '24px' } }}
            >
              <Form
                form={form}
                layout="vertical"
                onValuesChange={(changedValues, allValues) => {
                  // Jeśli zmieniła się tylko sekwencja, nie szukaj automatycznie
                  if (Object.keys(changedValues).length === 1 && 'sequence' in changedValues) {
                    return;
                  }
                  // Jeśli zmieniły się suwaki, czekamy na onAfterChange, aby nie przeciążać bazy
                  if ('angle_range' in changedValues || 'nt_range' in changedValues) {
                    return;
                  }
                  onSearch(allValues);
                }}
              >
                <Row gutter={[16, 0]}>
                  <Col xs={24} sm={12} md={6} lg={3}>
                    <Form.Item name="search_pdb" label="CIF ID">
                      <Select
                        placeholder="CIF ID"
                        allowClear
                        showSearch
                        filterOption={(input, option) =>
                          (option?.value ?? '').toString().toLowerCase().includes(input.toLowerCase())
                        }
                        options={availableIds.map(id => ({ value: id, label: id }))}
                      />
                    </Form.Item>
                  </Col>
                  <Col xs={24} sm={12} md={6} lg={4}>
                    <Form.Item
                      name="sequence"
                      label="Sequence Pattern"
                      getValueFromEvent={(e) => e.target.value.toUpperCase().replace(/[^ACGU]/g, '')}
                    >
                      <Input.Search 
                        placeholder="e.g. GAAA" 
                        allowClear 
                        onSearch={() => onSearch(form.getFieldsValue())}
                        enterButton
                      />
                    </Form.Item>
                  </Col>
                  <Col xs={24} sm={12} md={12} lg={4}>
                   <Form.Item name="segment_type" label="Type">
                     <Select
                       mode="multiple"
                       placeholder="Select types"
                       allowClear
                       maxTagCount="responsive"
                       popupRender={(menu) => (
                         <>
                           {menu}
                           <Divider style={{ margin: '8px 0' }} />
                           <Space style={{ padding: '0 8px 4px' }}>
                             <Button type="link" size="small" onClick={() => {
                               const vals = availableTypes.map(t => t.value);
			       form.setFieldsValue({ segment_type: vals });
                               onSearch(form.getFieldsValue());
                             }}>
                               Select All
                             </Button>
                             <Button type="link" size="small" onClick={() => {
                               form.setFieldsValue({ segment_type: [] });
                               onSearch(form.getFieldsValue());
                             }}>
                               Deselect All
                             </Button>
                           </Space>
                         </>
                       )}
                     >
		       {availableTypes.map(type => (
 			 <Select.Option key={type.value} value={type.value}>{type.label}</Select.Option>
		       ))}
                     </Select>
                   </Form.Item>
                  </Col>
                  <Col xs={24} sm={12} md={12} lg={4}>
                    <Form.Item label="Coaxially Stacked">
                      <Space.Compact style={{ width: '100%' }}>
                        <Form.Item name="stacking_stem1" noStyle>
                          <Select 
                            placeholder="Stem A" 
                            allowClear 
                            options={stemOptions} 
                            style={{ width: '50%' }} 
                          />
                        </Form.Item>
                        <Form.Item name="stacking_stem2" noStyle>
                          <Select 
                            placeholder="Stem B" 
                            allowClear 
                            options={stemOptions} 
                            style={{ width: '50%' }} 
                          />
                        </Form.Item>
                      </Space.Compact>
                    </Form.Item>
                  </Col>
                  <Col xs={12} sm={12} md={6} lg={3}>
                   <Form.Item name="angle_range" label="Angle (0-50°)" initialValue={[0, 50]}>
                     <Slider 
                       range 
                       min={0} 
                       max={50} 
                       step={0.01}
                       onAfterChange={() => onSearch(form.getFieldsValue())}
                     />
                   </Form.Item>
                  </Col>
		  <Col xs={12} sm={12} md={6} lg={3}>
                   <Form.Item name="nt_range" label={`NT Range`} initialValue={[0, maxNtLimit]}>
                     <Slider 
                       range 
                       min={0} 
                       max={maxNtLimit} 
                       onAfterChange={() => onSearch(form.getFieldsValue())}
                     />
                   </Form.Item>
                  </Col>
                  <Col xs={24} lg={3} style={{ display: 'flex', alignItems: 'flex-end', justifyContent: 'flex-end', paddingBottom: '24px' }}>
                   <Button icon={<ClearOutlined />} onClick={clearFilters} size="large">
                     Reset Filters
                   </Button>
                  </Col>
                  </Row>
                  </Form>            </Card>
          </Col>

          <Col span={24}>
            <Row gutter={[16, 16]}>
              <Col xs={24} sm={24} md={5}>
                <Card style={{ height: '250px', borderRadius: '8px', display: 'flex', alignItems: 'center', justifyContent: 'center' }}>
                  <Statistic 
                    title="Results" 
                    value={totalResults} 
                    valueStyle={{ fontSize: '28px', color: '#1890ff' }} 
                  />
                </Card>
              </Col>
              <Col xs={24} sm={24} md={12}>
                <Card title={<span style={{ fontSize: '13px', color: '#888' }}>Bend Angle Distribution</span>} style={{ height: '250px', borderRadius: '8px' }} styles={{ body: { padding: '10px 15px' } }}>
                  <div style={{ width: '100%', height: '190px', overflowX: 'auto', overflowY: 'hidden' }}>
                    <div style={{ minWidth: '600px', height: '180px' }}>
                      <ResponsiveContainer width="100%" height="100%">
                        <BarChart 
                          data={prepareHistogramData()}
                          margin={{ top: 5, right: 10, left: -20, bottom: 30 }}
                          onClick={(state) => {
                            if (state && state.activeLabel) {
                              const dataPoint = prepareHistogramData().find(d => d.range === state.activeLabel);
                              if (dataPoint) handleBarClick(dataPoint);
                            }
                          }}
                          style={{ cursor: 'pointer' }}
                        >
                          <CartesianGrid strokeDasharray="3 3" vertical={false} />
                          <XAxis dataKey="range" fontSize={10} tickLine={false} />
                          <YAxis fontSize={10} axisLine={false} tickLine={false} />
                          <ChartTooltip 
                            cursor={{fill: 'rgba(0, 0, 0, 0.05)'}} 
                            wrapperStyle={{ zIndex: 9999 }}
                            contentStyle={{ 
                              borderRadius: '8px', 
                              border: 'none', 
                              boxShadow: '0 4px 12px rgba(0,0,0,0.15)', 
                              fontSize: '11px'
                            }}
                            itemSorter={(item) => -(item.value as number)}
                            formatter={(value: number, name: string) => value > 0 ? [value, name] : [null, null]}
                            filterNull={true}
                          />
                          <Legend 
                            iconType="circle" 
                            wrapperStyle={{ 
                              fontSize: '10px',
                              paddingTop: '10px',
                              zIndex: 1,
                              pointerEvents: 'none'
                            }} 
                          />
                          {Object.entries(TYPE_COLORS).map(([type, color]) => (
                            <Bar 
                              key={type} 
                              dataKey={type} 
                              fill={color} 
                              fillOpacity={0.8}
                              barSize={15}
                            />
                          ))}
                        </BarChart>
                      </ResponsiveContainer>
                    </div>
                  </div>
                </Card>
              </Col>
              <Col xs={24} sm={24} md={7}>
                <Card title={<span style={{ fontSize: '13px', color: '#888' }}>Composition</span>} style={{ height: '250px', borderRadius: '8px' }} styles={{ body: { padding: '5px 15px' } }}>
                  <div style={{ width: '100%', height: '180px', minWidth: 0 }}>
                    <ResponsiveContainer width="100%" height="100%">
                      <PieChart>
                        <Pie
                          data={preparePieData()}
                          cx="50%"
                          cy="45%"
                          innerRadius={30}
                          outerRadius={50}
                          paddingAngle={5}
                          dataKey="value"
                        >
                          {preparePieData().map((entry, index) => (
                            <Cell 
                              key={`cell-${index}`} 
                              fill={TYPE_COLORS[entry.name] || PIE_COLORS[index % PIE_COLORS.length]} 
                            />
                          ))}
                        </Pie>
                        <ChartTooltip />
                        <Legend iconType="circle" wrapperStyle={{ fontSize: '10px', marginTop: '-10px' }} />
                      </PieChart>
                    </ResponsiveContainer>
                  </div>
                </Card>
              </Col>
            </Row>
          </Col>

          <Col span={24}>
            <Card style={{ borderRadius: '8px' }}>
              <div style={{ display: 'flex', justifyContent: 'space-between', alignItems: 'center', marginBottom: '16px', flexWrap: 'wrap', gap: '10px' }}>
                <Text strong style={{ fontSize: '16px' }}>
                  Found: <Tag color="blue">{totalResults}</Tag> results
                  <Text type="secondary" style={{ marginLeft: '15px', fontSize: '12px', fontWeight: 'normal' }}>
                    (Click any row to see details)
                  </Text>
                </Text>
                <Space wrap>
                  <Button 
                    icon={<FileTextOutlined />} 
                    onClick={exportToCsv}
                    disabled={results.length === 0}
                    type="default"
                  >
                    Export to CSV
                  </Button>
                  <Button 
                    icon={<DownloadOutlined />} 
                    onClick={handleBatchDownload}
                    disabled={results.length === 0}
                    type="primary"
                  >
                    Download All Results (.zip)
                  </Button>
                </Space>
              </div>
              <Table 
                dataSource={results} 
                columns={columns} 
                rowKey={(record) => `${record.type}-${record.id}`}
                loading={loading}
                pagination={{ 
                  current: currentPage,
                  pageSize: pageSize,
                  total: totalResults,
                  showSizeChanger: true,
                  pageSizeOptions: ['10', '20', '50', '100'],
                  position: ['bottomCenter']
                }}
                onChange={(pagination, filters, sorter: any) => {
                  onSearch(form.getFieldsValue(), pagination.current, pagination.pageSize, sorter);
                }}
                bordered
                size="middle"
                scroll={{ x: 800 }}
                expandable={{
                  expandedRowRender,
                  expandRowByClick: true,
                }}
                onRow={() => ({
                  style: { cursor: 'pointer' }
                })}
              />
            </Card>
          </Col>
        </Row>
      } />
      <Route path="/help" element={<HelpPage />} />
      <Route path="/about" element={<AboutPage />} />
      <Route path="/cite" element={<CitePage />} />
      <Route path="*" element={<Link to="/"><Button type="primary">Go Home</Button></Link>} />
    </Routes>
  </Content>
      <Footer style={{ textAlign: 'center', color: '#888', background: '#fff', borderTop: '1px solid #e8e8e8' }}>
        <strong>RNAbridge</strong> ©2026 | System for local RNA geometry analysis
      </Footer>

      <style>{`
        .preview-modal-container {
          display: flex;
          flex-direction: row;
          gap: 15px;
          min-height: 730px;
        }
        .preview-view-3d, .preview-view-2d {
          flex: 1;
          background: #ffffff;
          border-radius: 8px;
          border: 1px solid #f0f0f0;
          padding: 15px;
          display: flex;
          flex-direction: column;
        }
        @media (max-width: 1200px) {
          .preview-modal-container {
            display: none;
          }
          .preview-modal-tabs {
            display: block;
          }
        }
        @media (min-width: 1201px) {
          .preview-modal-tabs {
            display: none;
          }
        }
        `}</style>




        <Modal
        title={
          <Space>
            <DesktopOutlined />
            <Text strong>Structure Preview: {currentRecord?.pdb_id}</Text>
            <Tag color="blue">{currentRecord?.segment_count_folder}</Tag>
          </Space>
        }
        open={combinedModalVisible}
        onCancel={() => setCombinedModalVisible(false)}
        footer={[
          <Button key="close" type="primary" onClick={() => setCombinedModalVisible(false)}>Close</Button>
        ]}
        width="90%"
        style={{ top: 20, maxWidth: '1500px' }}
        centered={false}
        destroyOnClose
        >


        {currentRecord && (
          <>
            {/* LARGE SCREENS: Side-by-Side */}
            <div className="preview-modal-container">
              <div className="preview-view-2d">
                <Title level={5} style={{ textAlign: 'center', marginBottom: '15px' }}><EyeOutlined /> 2D Schematic</Title>
                <div style={{ height: '680px', width: '100%', display: 'flex', alignItems: 'center', justifyContent: 'center', overflow: 'hidden' }}>
                  {currentRecord.path_svg ? (
                    <Image 
                      src={`${API_BASE}/api/files/${currentRecord.path_svg}`} 
                      style={{ width: '100%', height: '680px', objectFit: 'contain' }} 
                      preview={false}
                    />
                  ) : <div style={{ color: '#888' }}>No 2D visualization available</div>}
                </div>
              </div>
              <div className="preview-view-3d">
                <Title level={5} style={{ textAlign: 'center', marginBottom: '15px' }}><DeploymentUnitOutlined /> 3D Model</Title>
                <Suspense fallback={<div style={{ textAlign: 'center', padding: '100px' }}><Spin size="large">Loading 3D model...</Spin></div>}>
                  <MolstarViewer url={`${API_BASE}/api/files/${currentRecord.path_cif}`} highlights={highlightData} height="680px" />
                </Suspense>
              </div>
            </div>

            {/* SMALL SCREENS: Tabs */}
            <div className="preview-modal-tabs">
              <Tabs
                defaultActiveKey="2d"
                centered
                items={[
                  {
                    key: '2d',
                    label: <span><EyeOutlined /> 2D Schematic</span>,
                    children: (
                      <div style={{ height: '700px', width: '100%', display: 'flex', alignItems: 'center', justifyContent: 'center', overflow: 'hidden', background: '#fff', borderRadius: '8px', border: '1px solid #f0f0f0' }}>
                        {currentRecord.path_svg ? (
                          <img 
                            src={`${API_BASE}/api/files/${currentRecord.path_svg}`} 
                            style={{ width: '100%', height: '100%', objectFit: 'contain' }} 
                            alt="2D Schematic"
                          />
                        ) : <div style={{ color: '#888' }}>No 2D visualization available</div>}
                      </div>
                    ),
                  },
                  {
                    key: '3d',
                    label: <span><DeploymentUnitOutlined /> 3D Model</span>,
                    children: (
                      <div style={{ height: '700px', background: '#fff', borderRadius: '8px', border: '1px solid #f0f0f0', overflow: 'hidden' }}>
                        <Suspense fallback={<div style={{ textAlign: 'center', padding: '100px' }}><Spin size="large">Loading 3D model...</Spin></div>}>
                          <MolstarViewer url={`${API_BASE}/api/files/${currentRecord.path_cif}`} highlights={highlightData} height="680px" />
                        </Suspense>
                      </div>
                    ),
                  },
                ]}
              />
            </div>

            {/* LEGEND SECTION */}
            <div style={{ marginTop: '20px', padding: '15px', background: '#f9f9f9', borderRadius: '8px', border: '1px solid #eee' }}>
              <Row gutter={[24, 16]}>
                <Col xs={24} md={8}>
                  <Text strong style={{ display: 'block', marginBottom: '8px' }}>Line Notations</Text>
                  <Space direction="vertical" size={4}>
                    <Space><div style={{ width: '20px', height: '2px', background: '#1890ff' }}></div> <Text style={{ fontSize: '12px' }}>Non-canonical base pairs</Text></Space>
                    <Space><div style={{ width: '20px', height: '2px', background: '#ff0000' }}></div> <Text style={{ fontSize: '12px' }}>Stacking paths</Text></Space>
                  </Space>
                </Col>
                <Col xs={24} md={8}>
                  <Text strong style={{ display: 'block', marginBottom: '8px' }}>Stem Highlights</Text>
                  <Space direction="vertical" size={4}>
                    <Space><div style={{ width: '12px', height: '12px', background: '#000000', borderRadius: '2px' }}></div> <Text style={{ fontSize: '12px' }}>Coaxial Stacking</Text></Space>
                    <Space><div style={{ width: '12px', height: '12px', background: '#d3d3d3', borderRadius: '2px' }}></div> <Text style={{ fontSize: '12px' }}>Generic / Other Stems</Text></Space>
                  </Space>
                </Col>
                <Col xs={24} md={8}>
                  <Text strong style={{ display: 'block', marginBottom: '8px' }}>Component Palettes</Text>
                  <Space direction="vertical" size={8}>
                    <div>
                      <Text type="secondary" style={{ display: 'block', marginBottom: '4px', fontSize: '12px' }}>Motifs:</Text>
                      <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px' }}>
                        {["#1890FF", "#A0D911", "#722ED1", "#13C2C2", "#EB2F96", "#237804", "#FAAD14", "#2F54EB"].map(c => (
                          <div key={c} style={{ width: '12px', height: '12px', background: c, borderRadius: '2px' }}></div>
                        ))}
                      </div>
                    </div>
                    <div>
                      <Text type="secondary" style={{ display: 'block', marginBottom: '4px', fontSize: '12px' }}>Unpaired Strands Linking Junction Stems:</Text>
                      <div style={{ display: 'flex', flexWrap: 'wrap', gap: '4px' }}>
                        {["#FA8C16", "#FADB14", "#EB2F96", "#FA541C", "#A0D911", "#D4380D", "#C41D7F", "#722ED1"].map(c => (
                          <div key={c} style={{ width: '12px', height: '12px', background: c, borderRadius: '2px' }}></div>
                        ))}
                      </div>
                    </div>
                  </Space>
                </Col>
              </Row>
            </div>
          </>
        )}
      </Modal>

      <Modal
        title={<span><FileTextOutlined /> Raw Data (JSON)</span>}
        open={jsonPreviewVisible}
        onCancel={() => setJsonPreviewVisible(false)}
        footer={[
          <Button key="close" type="primary" onClick={() => setJsonPreviewVisible(false)}>Close</Button>
        ]}
        width={800}
        centered
      >
        <div style={{ background: '#1e1e1e', color: '#d4d4d4', padding: '15px', borderRadius: '4px', maxHeight: '500px', overflowY: 'auto' }}>
          <pre style={{ margin: 0, fontSize: '12px', fontFamily: 'monospace' }}>{currentJsonContent}</pre>
        </div>
      </Modal>

      <Modal
        title={<span><FileTextOutlined /> CSV File Preview</span>}
        open={csvPreviewVisible}
        onCancel={() => setCsvPreviewVisible(false)}
        width={1000}
        footer={[
          <Button key="cancel" onClick={() => setCsvPreviewVisible(false)}>
            Cancel
          </Button>,
          <Button key="download" type="primary" icon={<DownloadOutlined />} onClick={downloadCsv}>
            Download CSV file
          </Button>
        ]}
      >
        <div style={{ background: '#f5f5f5', padding: '15px', borderRadius: '4px', maxHeight: '600px', overflowY: 'auto' }}>
          <pre style={{ margin: 0, fontSize: '13px', whiteSpace: 'pre-wrap', fontFamily: 'monospace' }}>{csvContent}</pre>
        </div>
        <Text type="secondary" style={{ display: 'block', marginTop: '10px' }}>
          The file contains {results.length} rows of data (plus header).
        </Text>
      </Modal>
    </Layout>
  );
};

export default App;
