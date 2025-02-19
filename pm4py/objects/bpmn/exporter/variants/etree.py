'''
    This file is part of PM4Py (More Info: https://pm4py.fit.fraunhofer.de).

    PM4Py is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    PM4Py is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with PM4Py.  If not, see <https://www.gnu.org/licenses/>.
'''
import uuid

from pm4py.objects.bpmn.obj import BPMN
from pm4py.util import constants, exec_utils
from enum import Enum


class Parameters(Enum):
    ENCODING = "encoding"


def apply(bpmn_graph, target_path, parameters=None):
    """
    Exports the BPMN diagram to a file

    Parameters
    -------------
    bpmn_graph
        BPMN diagram
    target_path
        Target path
    parameters
        Possible parameters of the algorithm
    """
    xml_string = get_xml_string(bpmn_graph, parameters=parameters)
    F = open(target_path, "wb")
    F.write(xml_string)
    F.close()


def get_xml_string(bpmn_graph, parameters=None):
    if parameters is None:
        parameters = {}

    encoding = exec_utils.get_param_value(Parameters.ENCODING, parameters, constants.DEFAULT_ENCODING)

    layout = bpmn_graph.get_layout()

    import xml.etree.ElementTree as ET
    from xml.dom import minidom

    definitions = ET.Element("bpmn:definitions")
    definitions.set("xmlns:bpmn", "http://www.omg.org/spec/BPMN/20100524/MODEL")
    definitions.set("xmlns:bpmndi", "http://www.omg.org/spec/BPMN/20100524/DI")
    definitions.set("xmlns:omgdc", "http://www.omg.org/spec/DD/20100524/DC")
    definitions.set("xmlns:omgdi", "http://www.omg.org/spec/DD/20100524/DI")
    definitions.set("xmlns:xsi", "http://www.w3.org/2001/XMLSchema-instance")
    definitions.set("targetNamespace", "http://www.signavio.com/bpmn20")
    definitions.set("typeLanguage", "http://www.w3.org/2001/XMLSchema")
    definitions.set("expressionLanguage", "http://www.w3.org/1999/XPath")
    definitions.set("xmlns:xsd", "http://www.w3.org/2001/XMLSchema")

    all_processes = set()
    process_planes = {}
    process_process = {}
    for node in bpmn_graph.get_nodes():
        all_processes.add(node.get_process())

    for process in all_processes:
        p = ET.SubElement(definitions, "bpmn:process",
                          {"id": "id" + process, "isClosed": "false", "isExecutable": "false",
                           "processType": "None"})
        process_process[process] = p

    diagram = ET.SubElement(definitions, "bpmndi:BPMNDiagram", {"id": "id" + str(uuid.uuid4()), "name": "diagram"})

    for process in all_processes:
        plane = ET.SubElement(diagram, "bpmndi:BPMNPlane",
                              {"bpmnElement": "id" + process, "id": "id" + str(uuid.uuid4())})
        process_planes[process] = plane

    for node in bpmn_graph.get_nodes():
        process = node.get_process()

        node_shape = ET.SubElement(process_planes[process], "bpmndi:BPMNShape",
                                   {"bpmnElement": node.get_id(), "id": node.get_id() + "_gui"})
        node_shape_layout = ET.SubElement(node_shape, "omgdc:Bounds",
                                          {"height": str(layout.get(node).get_height()), "width": str(layout.get(node).get_width()),
                                           "x": str(layout.get(node).get_x()),
                                           "y": str(layout.get(node).get_y())})

    for flow in bpmn_graph.get_flows():
        process = flow.get_process()

        flow_shape = ET.SubElement(process_planes[process], "bpmndi:BPMNEdge",
                                   {"bpmnElement": "id" + str(flow.get_id()),
                                    "id": "id" + str(flow.get_id()) + "_gui"})
        for x, y in layout.get(flow).get_waypoints():
            waypoint = ET.SubElement(flow_shape, "omgdi:waypoint", {"x": str(x), "y": str(y)})

    for node in bpmn_graph.get_nodes():
        process = process_process[node.get_process()]

        if isinstance(node, BPMN.StartEvent):
            isInterrupting = "true" if node.get_isInterrupting() else "false"
            parallelMultiple = "true" if node.get_parallelMultiple() else "false"
            task = ET.SubElement(process, "bpmn:startEvent",
                                 {"id": node.get_id(), "isInterrupting": isInterrupting, "name": node.get_name(),
                                  "parallelMultiple": parallelMultiple})
        elif isinstance(node, BPMN.EndEvent):
            task = ET.SubElement(process, "bpmn:endEvent", {"id": node.get_id(), "name": node.get_name()})
        elif isinstance(node, BPMN.IntermediateCatchEvent):
            task = ET.SubElement(process, "bpmn:intermediateCatchEvent", {"id": node.get_id(), "name": node.get_name()})
        elif isinstance(node, BPMN.IntermediateThrowEvent):
            task = ET.SubElement(process, "bpmn:intermediateThrowEvent", {"id": node.get_id(), "name": node.get_name()})
        elif isinstance(node, BPMN.BoundaryEvent):
            task = ET.SubElement(process, "bpmn:boundaryEvent", {"id": node.get_id(), "name": node.get_name()})
        elif isinstance(node, BPMN.Task):
            task = ET.SubElement(process, "bpmn:task", {"id": node.get_id(), "name": node.get_name()})
        elif isinstance(node, BPMN.SubProcess):
            task = ET.SubElement(process, "bpmn:subProcess", {"id": node.get_id(), "name": node.get_name()})
        elif isinstance(node, BPMN.ExclusiveGateway):
            task = ET.SubElement(process, "bpmn:exclusiveGateway",
                                 {"id": node.get_id(), "gatewayDirection": node.get_gateway_direction().value,
                                  "name": ""})
        elif isinstance(node, BPMN.ParallelGateway):
            task = ET.SubElement(process, "bpmn:parallelGateway",
                                 {"id": node.get_id(), "gatewayDirection": node.get_gateway_direction().value,
                                  "name": ""})
        elif isinstance(node, BPMN.InclusiveGateway):
            task = ET.SubElement(process, "bpmn:inclusiveGateway",
                                 {"id": node.get_id(), "gatewayDirection": node.get_gateway_direction().value,
                                  "name": ""})
        else:
            raise Exception("Unexpected node type.")

        for in_arc in node.get_in_arcs():
            arc_xml = ET.SubElement(task, "bpmn:incoming")
            arc_xml.text = "id" + str(in_arc.get_id())

        for out_arc in node.get_out_arcs():
            arc_xml = ET.SubElement(task, "bpmn:outgoing")
            arc_xml.text = "id" + str(out_arc.get_id())

    for flow in bpmn_graph.get_flows():
        process = process_process[flow.get_process()]

        source = flow.get_source()
        target = flow.get_target()
        flow_xml = ET.SubElement(process, "bpmn:sequenceFlow", {"id": "id" + str(flow.get_id()), "name": flow.get_name(),
                                                           "sourceRef": str(source.get_id()),
                                                           "targetRef": str(target.get_id())})

    return minidom.parseString(ET.tostring(definitions)).toprettyxml(encoding=encoding)
